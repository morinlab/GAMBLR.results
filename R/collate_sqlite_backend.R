#' @title Connect to the shared collated-results SQLite database.
#'
#' @description Returns a cached DBI connection to the SQLite database that
#' backs [GAMBLR.results::collate_results_db], creating the file (and its
#' parent directory) on first connection if it doesn't exist yet. Unlike
#' `gambl_mutations_db()` in `GAMBLR.data`, this database is writable and
#' local to this repo -- it isn't distributed via GitHub Releases, since
#' its contents (collated results) are internal/pre-publication.
#'
#' @param db_path Optional explicit path to the .db file. Defaults to
#' `results_merged$collated_db` under `project_base`, resolved via
#' `config.yml`.
#'
#' @return A DBIConnection.
#'
#' @import DBI RSQLite GAMBLR.helpers
#' @keywords internal
#' @noRd
gambl_collated_db <- function(db_path = NULL) {
  if (is.null(db_path)) {
    base <- GAMBLR.helpers::check_config_value(config::get("project_base"))
    rel_path <- GAMBLR.helpers::check_config_value(config::get("results_merged")$collated_db)
    db_path <- paste0(base, rel_path)
  }
  cached <- getOption("gamblr.collated.con")
  if (!is.null(cached) && DBI::dbIsValid(cached) &&
      identical(normalizePath(cached@dbname), normalizePath(db_path, mustWork = FALSE))) {
    return(cached)
  }
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  options(gamblr.collated.con = con)
  con
}

#' @title Map an R column's class to a SQLite column type affinity.
#'
#' @description Used only when a column needs to be added to an existing
#' table via `ALTER TABLE ... ADD COLUMN`, since that requires an explicit
#' type (unlike `DBI::dbWriteTable()`'s own type inference on initial table
#' creation).
#'
#' @param x A vector (one column's worth of data).
#'
#' @return A single SQLite type affinity string.
#'
#' @keywords internal
#' @noRd
sqlite_type_affinity <- function(x) {
  if (is.integer(x) || is.logical(x)) return("INTEGER")
  if (is.numeric(x)) return("REAL")
  "TEXT"
}

#' @title Get the set of keys already present in a collate results table.
#'
#' @description Returns the distinct `key_cols` combinations already stored
#' in `table_name`, or a zero-row placeholder if the table doesn't exist
#' yet (i.e. nothing is "already present" -- everything is missing).
#'
#' @param con A DBIConnection, from `gambl_collated_db()`.
#' @param table_name Name of the table to check.
#' @param key_cols Columns identifying a unique result. Default
#' `c("sample_id", "seq_type")`.
#'
#' @return A data frame with just `key_cols`.
#'
#' @import DBI dplyr
#' @keywords internal
#' @noRd
get_existing_collate_keys <- function(con, table_name, key_cols = c("sample_id", "seq_type")) {
  if (!DBI::dbExistsTable(con, table_name)) {
    empty <- as.data.frame(matrix(character(0), ncol = length(key_cols)))
    names(empty) <- key_cols
    return(empty)
  }
  cols_sql <- paste(DBI::dbQuoteIdentifier(con, key_cols), collapse = ", ")
  query <- paste0("SELECT DISTINCT ", cols_sql, " FROM ", DBI::dbQuoteIdentifier(con, table_name))
  DBI::dbGetQuery(con, query)
}

#' @title Write (create or upsert into) a collate results table.
#'
#' @description Creates `table_name` if it doesn't exist yet (schema
#' inferred from `new_data`). If it does exist, adds any columns present in
#' `new_data` but missing from the table (`ALTER TABLE ... ADD COLUMN`,
#' since a `collate_*` helper's output can gain columns over time), then
#' replaces any existing rows sharing a `key_cols` combination with
#' `new_data` and inserts the rest -- delete-matching-rows-then-insert,
#' rather than a raw `INSERT ... ON CONFLICT`, so this doesn't depend on a
#' unique index existing on `key_cols`. Wrapped in a transaction so a
#' failure partway through leaves the table in its prior state rather than
#' a half-written one.
#'
#' @param con A DBIConnection, from `gambl_collated_db()`.
#' @param table_name Name of the table to write to.
#' @param new_data A data frame including `key_cols` plus whatever result
#' columns this `collate_*` helper produced. Must contain at least one row.
#' @param key_cols Columns identifying a unique result. Default
#' `c("sample_id", "seq_type")`.
#'
#' @return Invisibly, `new_data`.
#'
#' @import DBI dplyr
#' @keywords internal
#' @noRd
write_collate_table <- function(con, table_name, new_data, key_cols = c("sample_id", "seq_type")) {
  stopifnot(nrow(new_data) > 0, all(key_cols %in% names(new_data)))

  if (!DBI::dbExistsTable(con, table_name)) {
    DBI::dbWriteTable(con, table_name, new_data)
    return(invisible(new_data))
  }

  existing_cols <- DBI::dbListFields(con, table_name)
  new_cols <- setdiff(names(new_data), existing_cols)
  for (col in new_cols) {
    ddl <- paste0(
      "ALTER TABLE ", DBI::dbQuoteIdentifier(con, table_name),
      " ADD COLUMN ", DBI::dbQuoteIdentifier(con, col),
      " ", sqlite_type_affinity(new_data[[col]])
    )
    DBI::dbExecute(con, ddl)
  }

  DBI::dbBegin(con)
  tryCatch({
    key_combos <- unique(new_data[key_cols])
    # A single DELETE whose WHERE clause OR's together one condition per key
    # combo hits SQLite's default max expression tree depth (1000) once
    # new_data covers more than a few hundred keys -- e.g. a registry entry
    # marked batchable = FALSE (see collate_registry.R) writes its entire
    # "missing" set in one write_collate_table() call, which can easily be
    # thousands of rows. Chunk the DELETE well under that limit; each chunk
    # still runs inside the same transaction, so this stays atomic.
    delete_chunk_size <- 200
    chunk_starts <- seq(1, nrow(key_combos), by = delete_chunk_size)
    for (start in chunk_starts) {
      end <- min(start + delete_chunk_size - 1, nrow(key_combos))
      chunk <- key_combos[start:end, , drop = FALSE]
      where_clauses <- apply(chunk, 1, function(row) {
        conditions <- mapply(function(col, val) {
          paste0(DBI::dbQuoteIdentifier(con, col), " = ", DBI::dbQuoteString(con, as.character(val)))
        }, key_cols, row)
        paste0("(", paste(conditions, collapse = " AND "), ")")
      })
      delete_sql <- paste0(
        "DELETE FROM ", DBI::dbQuoteIdentifier(con, table_name),
        " WHERE ", paste(where_clauses, collapse = " OR ")
      )
      DBI::dbExecute(con, delete_sql)
    }
    DBI::dbWriteTable(con, table_name, new_data, append = TRUE)
    DBI::dbCommit(con)
  }, error = function(e) {
    DBI::dbRollback(con)
    stop(e)
  })

  invisible(new_data)
}

# --- provenance -------------------------------------------------------------------------------------------------

# Details of each input file (path, modification time in UTC, size, md5); NA for files that don't exist
collate_input_details <- function(inputs, md5 = TRUE) {
  exists <- !is.na(inputs$input_path) & file.exists(inputs$input_path)
  info <- file.info(ifelse(exists, inputs$input_path, NA_character_))
  data.frame(
    sample_id = inputs$sample_id,
    input_path = inputs$input_path,
    input_mtime = ifelse(exists, format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), NA_character_),
    input_size = ifelse(exists, as.numeric(info$size), NA_real_),
    input_md5 = if (md5) ifelse(exists, unname(tools::md5sum(ifelse(exists, inputs$input_path, ""))), NA_character_)
                else NA_character_,
    stringsAsFactors = FALSE
  )
}

# The settings a registry entry's results depend on, as one-row character values (empty list if none)
collate_provenance_values <- function(entry) {
  if (is.null(entry$provenance)) return(list())
  values <- entry$provenance()
  lapply(values, function(v) paste(as.character(v), collapse = ","))
}

# Add provenance to rows about to be written: when (computed_at, UTC) and with what (gamblr_results_version and the
# entry's settings), and, for entries with input_files, the input file's path, modification time, size and md5.
stamp_collate_provenance <- function(new_cols, entry, scope) {
  new_cols$computed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  new_cols$gamblr_results_version <- as.character(utils::packageVersion("GAMBLR.results"))
  for (nm in names(prov <- collate_provenance_values(entry))) new_cols[[nm]] <- prov[[nm]]
  if (!is.null(entry$input_files)) {
    details <- collate_input_details(do.call(entry$input_files, list(scope)))
    new_cols <- dplyr::left_join(dplyr::select(new_cols, -dplyr::any_of(setdiff(names(details), "sample_id"))),
                                 details, by = "sample_id")
  }
  new_cols
}

# Sample_ids of `scope` whose stored rows in `table_name` are stale: written before provenance was recorded (no
# computed_at), computed with different settings than the entry's current ones, or, for entries with input_files,
# whose input file's modification time or size has changed (or is now missing). Cheap: no md5 is computed.
collate_stale_ids <- function(con, table_name, entry, scope, seq) {
  if (!DBI::dbExistsTable(con, table_name)) return(character(0))
  stored <- DBI::dbReadTable(con, table_name)
  stored <- stored[stored$seq_type == seq & stored$sample_id %in% scope$sample_id, , drop = FALSE]
  if (nrow(stored) == 0) return(character(0))
  stale <- if ("computed_at" %in% names(stored)) is.na(stored$computed_at) else rep(TRUE, nrow(stored))
  for (nm in names(prov <- collate_provenance_values(entry))) {
    stale <- stale | if (nm %in% names(stored)) is.na(stored[[nm]]) | stored[[nm]] != prov[[nm]] else TRUE
  }
  # only rows not already stale need their input file checked (rows written before provenance existed are stale
  # without looking at any file)
  check <- which(!stale)
  if (!is.null(entry$input_files) && length(check) > 0) {
    scope_rows <- scope[match(stored$sample_id[check], scope$sample_id), , drop = FALSE]
    now <- collate_input_details(do.call(entry$input_files, list(scope_rows)), md5 = FALSE)
    was_mtime <- stored$input_mtime[check]
    was_size <- as.numeric(stored$input_size[check])
    # a file that was missing when the row was computed and still is doesn't make the row stale; one that has
    # appeared, disappeared or changed does
    missing_now <- is.na(now$input_mtime)
    missing_then <- is.na(was_mtime)
    changed <- missing_now != missing_then
    both <- !missing_now & !missing_then
    changed[both] <- now$input_mtime[both] != was_mtime[both] | is.na(was_size[both]) |
      now$input_size[both] != was_size[both]
    stale[check] <- changed
  }
  unique(stored$sample_id[stale])
}

# Names of the provenance columns a registry entry's table can hold
collate_provenance_columns <- function(entry) {
  c("computed_at", "gamblr_results_version", names(collate_provenance_values(entry)),
    if (!is.null(entry$input_files)) c("input_path", "input_mtime", "input_size", "input_md5"))
}
