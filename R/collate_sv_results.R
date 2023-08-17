#' @title Collate SV Results.
#'
#' @description Determine and summarize which cases have specific oncogene SVs.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param tool Name of tool (optional, default is manta).
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param oncogenes Which oncogenes to collate SVs from.
#'
#' @return Data frame with additional columns ({tool}_{oncogene} and {tool}_{oncogene}_{partner}).
#'
#' @import dplyr
#'
#' @noRd
#'
#' @examples
#' results = collate_samples_sv_results(sample_table = samples,
#'                                      tool = "manta",
#'                                      oncogenes = c("MYC", "BCL2"))
#'
collate_sv_results = function(sample_table,
                              tool = "manta",
                              seq_type_filter = "genome",
                              oncogenes = c("MYC", "BCL2", "BCL6", "CCND1", "IRF4")){
  if(seq_type_filter!="genome"){
    message("skipping sv for this seq_type")
    return(sample_table)
  }
  if(tool == "manta"){
    all_svs = get_manta_sv()
  }
  annotated_svs = annotate_sv(all_svs) %>%
  dplyr::filter(!is.na(partner))

  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
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
  out_table = sample_table
  for(oncogene in oncogenes){
    out_table = multiout(out_table, annotated_svs, "manta", oncogene)
  }
  return(out_table)
}
