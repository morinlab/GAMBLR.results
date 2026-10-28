#' @title Compute oncogene SV summaries for a sample scope.
#'
#' @description Core computation behind `collate_sv_results()`, extracted
#' so it can also be called directly by [collate_results_db()] for just the
#' subset of samples missing from its cache table. SV calling only applies
#' to genome seq_type -- a batch with no genome samples in it (e.g. a
#' capture batch, since `collate_results_db()` loops over every seq_type
#' present in the caller's metadata) gets a full NA row per sample instead
#' of being skipped outright, so it's cached as "not applicable" rather
#' than re-attempted (reloading all SV calls) on every future run. This
#' inference is based on `sample_table$seq_type` rather than a static
#' `seq_type_filter` argument because the dispatcher's registered
#' `extra_args` are fixed per function, not re-evaluated per seq_type --
#' the data itself is the only thing that reliably varies per call.
#'
#' @param sample_table A data frame with `sample_id` (and, ideally,
#' `seq_type`) columns, scoping which samples to compute for.
#' @param tool Name of tool (optional, default is manta).
#' @param oncogenes Which oncogenes to collate SVs from.
#'
#' @return A data frame with `sample_id` plus `manta_{oncogene}_sv` /
#' `manta_{oncogene}_partner` for each oncogene (column names always use
#' the `manta_` prefix regardless of `tool`, matching the original
#' implementation).
#'
#' @import dplyr GAMBLR.utils
#'
#' @keywords internal
#' @noRd
compute_sv_results_core <- function(sample_table,
                                    tool = "svar",
                                    oncogenes = c("MYC", "BCL2", "BCL6", "CCND1", "IRF4")){

  expected_cols <- unlist(lapply(oncogenes, function(og) paste0("manta_", og, c("_sv", "_partner"))))

  if(!is.null(sample_table$seq_type) && !("genome" %in% sample_table$seq_type)){
    result <- dplyr::select(sample_table, sample_id)
    for(col in expected_cols) result[[col]] <- NA
    return(result)
  }

  if(tool == "manta"){
    all_svs = get_manta_sv()
  }else if(tool == "svar"){
    all_svs = get_combined_sv()
  }
  annotated_svs = GAMBLR.utils::annotate_sv(all_svs) %>%
  dplyr::filter(!is.na(partner))

  multiout = function(df,
                       annotated,
                       tool,
                       oncogene_name){

    some_fusions = dplyr::filter(annotated, gene == all_of(oncogene_name)) %>%
      group_by(tumour_sample_id) %>%
      arrange(partner) %>%
      dplyr::filter(row_number() == 1)

    df = mutate(df, "{tool}_{oncogene_name}_sv" := case_when(sample_id %in% some_fusions$tumour_sample_id ~ "POS", TRUE ~ "NEG"))
    some_fusions = some_fusions %>%
      dplyr::select(tumour_sample_id, partner) %>%
      mutate("{tool}_{oncogene_name}_partner" := partner) %>%
      dplyr::select(-partner)

    df = left_join(df, some_fusions, by = c("sample_id" = "tumour_sample_id"))
    return(df)
  }
  out_table = dplyr::select(sample_table, sample_id)
  for(oncogene in oncogenes){
    out_table = multiout(out_table, annotated_svs, "manta", oncogene)
  }
  return(out_table)
}

#' @title Collate SV Results.
#'
#' @description Determine and summarize which cases have specific oncogene SVs.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param tool Name of tool (optional, default is manta).
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param oncogenes Which oncogenes to collate SVs from.
#'
#' @return Data frame with additional columns ({tool}_{oncogene} and {tool}_{oncogene}_{partner}).
#'
#' @import dplyr GAMBLR.utils
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' results = collate_samples_sv_results(sample_table = samples,
#'                                      tool = "manta",
#'                                      oncogenes = c("MYC", "BCL2"))
#' }
collate_sv_results = function(sample_table,
                              tool = "svar",
                              seq_type_filter = "genome",
                              oncogenes = c("MYC", "BCL2", "BCL6", "CCND1", "IRF4")){
  if(seq_type_filter!="genome"){
    message("skipping sv for this seq_type")
    return(sample_table)
  }

  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>%
      dplyr::filter(seq_type=="genome") %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }

  new_cols <- compute_sv_results_core(
    sample_table = sample_table,
    tool = tool,
    oncogenes = oncogenes
  )
  out_table = left_join(sample_table, new_cols, by = "sample_id")

  return(out_table)
}
