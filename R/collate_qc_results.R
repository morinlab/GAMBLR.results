#' @title Collate Quality Control Results.
#'
#' @description Expand a metadata table horizontally with quality control metrics.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table df with sample ids in the first column.
#' @param seq_type_filter default is genome, capture is also available.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr readr glue GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#'   qc_metrics = collate_qc_results(sample_table = sample_table)
#' }
collate_qc_results = function(sample_table,
                              seq_type_filter = "genome"){

  if(! seq_type_filter %in% c("genome", "capture")){
    stop("Please provide a valid seq_type (\"genome\" or \"capture\").")
  }

  #get paths
  base = GAMBLR.helpers::check_config_value(config::get("project_base"))
  qc_template = GAMBLR.helpers::check_config_value(config::get("qc_met"))

  #icgc_dart
  unix_group = "icgc_dart"
  icgc_qc_path = glue::glue(qc_template)
  icgc_qc_path_full = paste0(base, icgc_qc_path)

  #gambl
  unix_group = "gambl"
  gambl_qc_path = glue::glue(qc_template)
  gambl_qc_path_full = paste0(base, gambl_qc_path)

  #read in icgc qc data, rename sample id column and filter on samples in sample id in sample_table
  icgc_qc = suppressMessages(read_tsv(icgc_qc_path_full)) %>%
      dplyr::rename(sample_id = UID) %>%
      dplyr::select(-SeqType)

  #read in gambl qc data, rename sample id column and filter on samples in sample id in sample_table
  gambl_qc = suppressMessages(read_tsv(gambl_qc_path_full)) %>%
    dplyr::rename(sample_id = UID) %>%
    dplyr::select(-SeqType)

  #join gambl and icgc QC data
  full_qc = rbind(gambl_qc, icgc_qc)
  sample_table = left_join(sample_table, full_qc)

  #print n samples with QC metrics
  qc_samples = length(unique(full_qc$sample_id))
  message(paste("QC metrics for", qc_samples, "samples retrieved."))

  return(sample_table)
}
