#' @title Collate GAMBL results via the SQLite-backed cache.
#'
#' @description The SQLite-backed counterpart to
#' [GAMBLR.results::collate_results]. For each function registered in the
#' internal `collate_registry`, computes results only for samples missing
#' from that function's table (or explicitly listed in `refresh`), then
#' returns one wide table joining every registered function's results onto
#' `these_samples_metadata`.
#'
#' @details This is additive, not a replacement: `collate_results()` and
#' the individual `collate_*_results()` wrapper functions are untouched and
#' keep working exactly as before, backed by the existing shared-TSV cache.
#' Only functions listed in the internal `collate_registry` are included
#' here -- see `R/collate_registry.R` for which ones and why.
#'
#' `these_samples_metadata` is deduplicated to one row per
#' `(sample_id, seq_type)` (keeping the first occurrence) before anything
#' else happens, with a message if it actually drops anything -- every
#' registered core function's output grain follows the metadata scope it's
#' given, so a duplicate row there would otherwise produce duplicate rows
#' in every table it touches, and in the final result.
#'
#' @param these_samples_metadata A metadata table with (at least)
#' `sample_id` and `seq_type` columns. Defaults to all genome and capture
#' samples from [GAMBLR.results::get_gambl_metadata] if omitted.
#' @param refresh Optional named list, keyed by registry name (e.g.
#' `"ssm_results"`), of sample_ids to force-recompute for that function
#' even if already present in its table. Functions not named here use the
#' default missing-only behaviour.
#' @param batch_size Compute at most this many samples per call to a
#' registered function's core function, writing each batch to its table
#' before moving to the next. Bounds memory to O(batch_size) rather than
#' O(samples missing) -- relevant since some core functions (e.g.
#' `compute_ssm_results_core()`, via `get_ssm_by_samples()`) hold every
#' batch member's data in memory at once before combining. Also means a
#' failure partway through a large catch-up run doesn't lose already-
#' completed batches: a re-run only recomputes what's left. Default 50.
#' Ignored for any function registered with `batchable = FALSE` in
#' `collate_registry.R` -- those compute their entire "missing" set in one
#' call regardless, since their cost comes from a small, fixed set of
#' shared file reads rather than scaling with sample count.
#' @param db_path Optional explicit path to the SQLite database, passed to
#' `gambl_collated_db()`.
#'
#' @return `these_samples_metadata` joined with every registered
#' function's result columns.
#'
#' @import dplyr DBI
#' @export
#'
#' @examples
#' \dontrun{
#' my_meta <- get_gambl_metadata() %>% dplyr::filter(pathology == "FL")
#' collated <- collate_results_db(these_samples_metadata = my_meta)
#'
#' # force ssm_results to recompute for two specific samples
#' collated <- collate_results_db(
#'   these_samples_metadata = my_meta,
#'   refresh = list(ssm_results = c("sample1", "sample2"))
#' )
#' }
collate_results_db <- function(these_samples_metadata, refresh = list(), batch_size = 50, db_path = NULL) {
  if (missing(these_samples_metadata)) {
    these_samples_metadata <- get_gambl_metadata() %>%
      dplyr::filter(seq_type %in% c("genome", "capture"))
  }
  if (!all(c("sample_id", "seq_type") %in% names(these_samples_metadata))) {
    stop("these_samples_metadata must include sample_id and seq_type columns.")
  }

  # Every core function keys its output on (sample_id, seq_type) via
  # dplyr::select(sample_table, sample_id) %>% left_join(...) -- if the
  # caller's own metadata has more than one row for the same
  # (sample_id, seq_type) (e.g. from an upstream join that fanned out),
  # that duplication flows straight through into what gets written to
  # each table, and then into the final joined result. Collapsing here,
  # once, keeps every downstream table and the final join clean
  # regardless of the cause. Not silent, since a duplicate row usually
  # means something upstream in the caller's metadata construction is
  # worth checking.
  n_before <- nrow(these_samples_metadata)
  these_samples_metadata <- dplyr::distinct(these_samples_metadata, sample_id, seq_type, .keep_all = TRUE)
  n_dropped <- n_before - nrow(these_samples_metadata)
  if (n_dropped > 0) {
    message(sprintf(
      "collate_results_db(): these_samples_metadata had %d duplicate (sample_id, seq_type) row(s); keeping the first occurrence of each. If unexpected, check how this metadata table was built.",
      n_dropped
    ))
  }

  con <- gambl_collated_db(db_path = db_path)
  seq_types <- unique(these_samples_metadata$seq_type)

  for (reg_name in names(collate_registry)) {
    entry <- collate_registry[[reg_name]]
    refresh_ids <- refresh[[reg_name]]
    if (is.null(refresh_ids)) refresh_ids <- character(0)

    for (seq in seq_types) {
      scope <- dplyr::filter(these_samples_metadata, seq_type == seq)
      requested_ids <- scope$sample_id

      existing <- get_existing_collate_keys(con, reg_name)
      existing_ids <- if (nrow(existing) > 0) {
        existing$sample_id[existing$seq_type == seq]
      } else {
        character(0)
      }

      to_compute <- union(
        setdiff(requested_ids, existing_ids),
        intersect(requested_ids, refresh_ids)
      )
      if (length(to_compute) == 0) next

      # Functions explicitly marked batchable = FALSE (see collate_registry.R)
      # have a cost dominated by a small, fixed set of shared file reads,
      # independent of how many samples are requested -- chunking those would
      # just re-read the same files once per batch for no benefit. Everything
      # else defaults to batchable, since most core functions' cost does
      # scale with the number of samples requested at once. A registered
      # entry may also override the batch size itself (entry$batch_size),
      # for a function whose per-sample cost warrants a smaller (or larger)
      # chunk than the global default.
      entry_batchable <- if (is.null(entry$batchable)) TRUE else isTRUE(entry$batchable)
      this_batch_size <- if (!entry_batchable) {
        length(to_compute)
      } else if (!is.null(entry$batch_size)) {
        entry$batch_size
      } else {
        batch_size
      }
      batches <- split(to_compute, ceiling(seq_along(to_compute) / this_batch_size))
      for (i in seq_along(batches)) {
        batch_ids <- batches[[i]]
        if (length(batches) > 1) {
          message(sprintf(
            "%s (%s): batch %d/%d (%d samples)",
            reg_name, seq, i, length(batches), length(batch_ids)
          ))
        }
        scope_subset <- dplyr::filter(scope, sample_id %in% batch_ids)
        call_args <- c(setNames(list(scope_subset), entry$metadata_arg), entry$extra_args)
        new_cols <- do.call(entry$core_fn, call_args)
        new_cols$seq_type <- seq
        write_collate_table(con, reg_name, new_cols)
      }
    }
  }

  result <- these_samples_metadata
  for (reg_name in names(collate_registry)) {
    if (!DBI::dbExistsTable(con, reg_name)) next
    table_data <- DBI::dbReadTable(con, reg_name)
    if (nrow(table_data) == 0) next

    # Skip any non-key column this table would produce that the caller's
    # own these_samples_metadata already has -- e.g. get_gambl_metadata()
    # already bakes QC columns into its output via its own internal
    # collate_qc_results() call (used for its min_corrected_cov filter),
    # so re-joining qc_results here would otherwise collide on identical
    # column names and get silently .x/.y suffixed by dplyr rather than
    # erroring. Whatever the caller's metadata already carries wins.
    value_cols <- setdiff(names(table_data), c("sample_id", "seq_type"))
    conflicting <- intersect(value_cols, names(result))
    if (length(conflicting) > 0) {
      message(sprintf(
        "collate_results_db(): %s column(s) already present in these_samples_metadata, skipping from the %s table: %s",
        length(conflicting), reg_name, paste(conflicting, collapse = ", ")
      ))
      table_data <- dplyr::select(table_data, -all_of(conflicting))
      value_cols <- setdiff(value_cols, conflicting)
    }
    if (length(value_cols) == 0) next

    result <- dplyr::left_join(result, table_data, by = c("sample_id", "seq_type"))
  }
  result
}
