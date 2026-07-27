#' @title Compute SSM summary statistics for a sample scope.
#'
#' @description Core computation behind `collate_ssm_results()`, extracted
#' so it can also be called directly by [collate_results_db()] for just
#' the subset of samples missing from its cache table. Takes a sample
#' scope in, returns only the new columns (`total_ssm`, `mean_vaf`,
#' `coding_ssm`), keyed by `sample_id` -- every `sample_id` present in
#' `sample_table` gets a row back, even if it had zero mutations (`NA`
#' counts), matching the join behaviour `collate_ssm_results()` always
#' had. Reads per-sample files via `get_ssm_by_samples()` rather than the
#' whole 200+MB merged MAF, since this is typically called for a small
#' "missing from the cache" scope, not the full cohort at once.
#'
#' @param sample_table A data frame with sample_id as a column, scoping
#' which samples to compute for.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param projection Specifies the projection, default is "grch37".
#' @param from_flatfile Accepted for signature compatibility with the
#' original implementation, but no longer branches on anything --
#' `from_flatfile = FALSE` never actually worked in the original (no
#' `muts` was ever defined for that case).
#' @param include_silent Logical parameter indicating whether to include
#' silent mutations into coding mutations. Default is FALSE.
#'
#' @return A data frame with `sample_id`, `total_ssm`, `mean_vaf`,
#' `coding_ssm`.
#'
#' @import dplyr
#'
#' @keywords internal
#' @noRd
compute_ssm_results_core <- function(sample_table,
                                     seq_type_filter = "genome",
                                     projection = "grch37",
                                     from_flatfile = TRUE,
                                     include_silent = FALSE){

  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }

  # Read per-sample files for just the requested scope (get_ssm_by_samples()
  # with its default subset_from_merge = FALSE), instead of loading the
  # entire merged MAF -- that file is 200+MB across the whole cohort, and
  # this function is now typically called for a handful of samples missing
  # from the cache, not the whole cohort at once (see collate_results_db()).
  # from_flatfile is accepted but no longer branches on anything: the
  # FALSE case never actually worked (no `muts` was ever defined for it).
  # min_read_support = 0 preserves this function's original behaviour of
  # not filtering by read support at all -- get_ssm_by_sample()'s own
  # filter is skipped when its threshold is 0/falsy.
  muts = get_ssm_by_samples(
    these_samples_metadata = sample_table,
    projection = projection,
    augmented = FALSE,
    min_read_support = 0,
    basic_columns = TRUE
  ) %>%
    dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,t_alt_count,t_ref_count)

  #get tally of total per sample
  muts = muts %>%
    dplyr::rename("sample_id" = "Tumor_Sample_Barcode") %>%
    # restrict to the requested scope before summarizing -- when called
    # for just a handful of missing samples (the common case once the
    # cache is mostly populated), this avoids tallying stats for
    # thousands of samples the caller didn't ask for.
    dplyr::filter(sample_id %in% sample_table$sample_id)

  muts = mutate(muts, vaf = t_alt_count/(t_alt_count + t_ref_count))
  muts_count = dplyr::select(muts, sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("total_ssm" = "n")

  muts_mean = muts %>%
    dplyr::select(sample_id, vaf) %>%
    group_by(sample_id) %>%
    summarize(mean_vaf = mean(vaf))

  coding_mut = dplyr::filter(muts, Variant_Classification %in% coding_class)
  coding_mut_count = coding_mut %>%
    dplyr::select(sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("coding_ssm" = "n")

  # every requested sample_id gets a row back, even with NA counts if it
  # had no mutations at all -- start from sample_table's own sample_ids,
  # not from whichever samples happened to show up in muts_count.
  result = dplyr::select(sample_table, sample_id) %>%
    left_join(muts_count, by = "sample_id") %>%
    left_join(muts_mean, by = "sample_id") %>%
    left_join(coding_mut_count, by = "sample_id")

  return(result)
}

#' @title Collate SSM Results.
#'
#' @description Compute summary statistics based on SSM calls.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param projection Specifies the projection, default is "grch37".
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations, default is TRUE.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr glue GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' ssm_results = colalte_ssm_results(sample_table = samples,
#'                                   include_silent = TRUE)
#' }
collate_ssm_results = function(sample_table,
                               seq_type_filter = "genome",
                               projection = "grch37",
                               from_flatfile = TRUE,
                               include_silent = FALSE){

  new_cols = compute_ssm_results_core(
    sample_table = sample_table,
    seq_type_filter = seq_type_filter,
    projection = projection,
    from_flatfile = from_flatfile,
    include_silent = include_silent
  )
  sample_table = dplyr::left_join(sample_table, new_cols, by = "sample_id")

  return(sample_table)
}
