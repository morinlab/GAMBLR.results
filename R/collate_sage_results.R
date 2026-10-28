#' @title Compute sage-flavour SSM summary statistics for a sample scope.
#'
#' @description Core computation behind `collate_sage_results()`, extracted
#' so it can also be called directly by [collate_results_db()] for just
#' the subset of samples missing from its cache table. Mirrors
#' `compute_ssm_results_core()` exactly, except it counts sage-flavour
#' calls instead of the default clustered/slms-3 ones (see `flavour` in
#' [GAMBLR.results::get_ssm_by_samples]). Uses `sage_*`-prefixed column
#' names deliberately: `ssm_results` and `sage_results` are both
#' registered in `collate_registry.R` and joined onto the same wide table
#' by `collate_results_db()`, and identically-named columns from two
#' registered tables would collide there (see the "already present,
#' skipping" behaviour in `collate_results_db()` -- the second table's
#' column would simply be dropped, not renamed).
#'
#' @param sample_table A data frame with sample_id as a column, scoping
#' which samples to compute for.
#' @param projection Specifies the projection, default is "grch37".
#' @param include_silent Logical parameter indicating whether to include
#' silent mutations into coding mutations. Default is FALSE.
#'
#' @return A data frame with `sample_id`, `sage_total_ssm`,
#' `sage_mean_vaf`, `sage_coding_ssm`.
#'
#' @import dplyr
#'
#' @keywords internal
#' @noRd
compute_sage_results_core <- function(sample_table,
                                      projection = "grch37",
                                      include_silent = FALSE){

  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }

  # Same rationale as compute_ssm_results_core(): read per-sample files for
  # just the requested scope rather than loading a whole merged MAF. sage
  # only has per-sample files anyway (subset_from_merge = TRUE isn't
  # supported for it), so this is the only option here, not just the
  # faster one.
  muts = get_ssm_by_samples(
    these_samples_metadata = sample_table,
    projection = projection,
    flavour = "sage",
    augmented = FALSE,
    min_read_support = 0,
    basic_columns = TRUE
  ) %>%
    dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,t_alt_count,t_ref_count)

  muts = muts %>%
    dplyr::rename("sample_id" = "Tumor_Sample_Barcode") %>%
    dplyr::filter(sample_id %in% sample_table$sample_id)

  muts = mutate(muts, vaf = t_alt_count/(t_alt_count + t_ref_count))
  muts_count = dplyr::select(muts, sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("sage_total_ssm" = "n")

  muts_mean = muts %>%
    dplyr::select(sample_id, vaf) %>%
    group_by(sample_id) %>%
    summarize(sage_mean_vaf = mean(vaf))

  coding_mut = dplyr::filter(muts, Variant_Classification %in% coding_class)
  coding_mut_count = coding_mut %>%
    dplyr::select(sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("sage_coding_ssm" = "n")

  # every requested sample_id gets a row back, even with NA counts if it
  # had no sage-flavour mutations at all.
  result = dplyr::select(sample_table, sample_id) %>%
    left_join(muts_count, by = "sample_id") %>%
    left_join(muts_mean, by = "sample_id") %>%
    left_join(coding_mut_count, by = "sample_id")

  return(result)
}

#' @title Collate sage-flavour SSM Results.
#'
#' @description Compute summary statistics based on sage-flavour SSM
#' calls, as a counterpart to `collate_ssm_results()` (which uses the
#' default clustered/slms-3 calls). See `flavour` in
#' [GAMBLR.results::get_ssm_by_samples] for what "sage" means.
#'
#' @details INTERNAL FUNCTION, primarily intended to be called via the
#' `sage_results` entry in `collate_registry.R` / [collate_results_db()],
#' not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param projection Specifies the projection, default is "grch37".
#' @param include_silent Logical parameter indicating whether to include
#' silent mutations into coding mutations. Default is FALSE.
#'
#' @return The sample table with additional sage_* columns.
#'
#' @import dplyr
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' sage_results = collate_sage_results(sample_table = samples)
#' }
collate_sage_results = function(sample_table,
                                projection = "grch37",
                                include_silent = FALSE){

  new_cols = compute_sage_results_core(
    sample_table = sample_table,
    projection = projection,
    include_silent = include_silent
  )
  sample_table = dplyr::left_join(sample_table, new_cols, by = "sample_id")

  return(sample_table)
}
