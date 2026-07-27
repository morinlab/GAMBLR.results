#' @title Compute mutational signature results for a sample scope.
#'
#' @description Core computation behind `collate_sbs_results()`, extracted
#' so it can also be called directly by [collate_results_db()] for just the
#' subset of samples missing from its cache table. SBS signature activities
#' were only computed for genome seq_type -- a batch with no genome samples
#' in it (e.g. a capture batch, since `collate_results_db()` loops over
#' every seq_type present in the caller's metadata) gets back just
#' `sample_id` with no signature columns at all, rather than being skipped
#' outright. That's still enough for `collate_results_db()` to cache those
#' sample rows as "not applicable" (any signature columns already on the
#' table are simply left `NULL` for them) instead of re-attempting the file
#' load on every future run. Unlike `compute_sv_results_core()`, the exact
#' set of signature column names isn't knowable without reading
#' `file_path`, so this skips that read entirely rather than constructing a
#' matching NA-filled row.
#'
#' @param sample_table A data frame with `sample_id` (and, ideally,
#' `seq_type`) columns, scoping which samples to compute for.
#' @param file_path Optional path to SBS file.
#' @param scale_vals Parameter not used?
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#'
#' @return A data frame with `sample_id` plus whatever signature columns
#' `file_path` contains (possibly none, for a non-genome scope).
#'
#' @import dplyr tibble GAMBLR.helpers
#'
#' @keywords internal
#' @noRd
compute_sbs_results_core <- function(sample_table,
                                     file_path,
                                     scale_vals = FALSE,
                                     sbs_manipulation = ""){

  if(!is.null(sample_table$seq_type) && !("genome" %in% sample_table$seq_type)){
    return(dplyr::select(sample_table, sample_id))
  }

  if(missing(file_path)){
    base = GAMBLR.helpers::check_config_value(config::get("project_base"))

    file_path = paste0(base,"icgc_dart/sigprofiler-1.0/02-extract/genome--hg38/BL_HGBL_DLBCL_FL_COMFL_CLL_MCL_B-ALL_PBL_DLBCL-BL-like_UNSPECIFIED_SCBC_MM_all/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt")
  }
  message(paste("loading",file_path))
  signatures = read.csv(file_path, sep = "\t", header = 1, row.names = 1)
  rs = rowSums(signatures)
  cn = colnames(signatures)
  new_sig = signatures
  if(sbs_manipulation == "scale"){
    sbs = apply(signatures, 2, function(x){x/rs}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
  }
  else if(sbs_manipulation == "log"){
    sbs = apply(signatures, 2, function(x){log(x + 1)}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
  }else{
    sbs = signatures %>%
      rownames_to_column("sample_id")
  }

  # every requested sample_id gets a row back, even one absent from the
  # signatures file (NA in the sig columns) -- left_join starts from the
  # full requested scope, not from whichever samples happen to show up in
  # `sbs`.
  result <- dplyr::select(sample_table, sample_id) %>%
    left_join(sbs, by = "sample_id")

  return(result)
}

#' @title Collate SBS Results.
#'
#' @description Bring in the results from mutational signature analysis.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param file_path Optional path to SBS file.
#' @param scale_vals Parameter not used?
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#'
#' @return A data frame with new columns added.
#'
#' @import dplyr tibble GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' collated = collate_sbs_results(sample_table = sample_table,
#'                                sbs_manipulation = sbs_manipulation)
#' }
collate_sbs_results = function(sample_table,
                               seq_type_filter = "genome",
                               file_path,
                               scale_vals = FALSE,
                               sbs_manipulation = ""){
  if(seq_type_filter!="genome"){
    message("skipping sbs for seq_type")
    return(sample_table)
  }

  new_cols <- compute_sbs_results_core(
    sample_table = sample_table,
    file_path = file_path,
    scale_vals = scale_vals,
    sbs_manipulation = sbs_manipulation
  )
  sample_table = left_join(sample_table, new_cols, by = "sample_id")

  return(sample_table)
}
