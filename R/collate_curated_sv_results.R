#' @title Collate Curated SV Results.
#'
#' @description Collate all SV calls from the genome data and summarize for main oncogenes of interest per sample.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return The sample table with additional columns.
#'
#' @import readr dplyr GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' gambl_results_derived = collate_curated_sv_results(gambl_results_derived)
#' }
collate_curated_sv_results = function(sample_table,
                                      seq_type_filter = "genome"){

  path_to_files = GAMBLR.helpers::check_config_value(config::get("derived_and_curated"))
  project_base = GAMBLR.helpers::check_config_value(config::get("project_base"))
  manual_files = dir(paste0(project_base, path_to_files), pattern = ".tsv")
  for(f in manual_files){
    full = paste0(project_base, path_to_files, f)
    this_data = suppressMessages(read_tsv(full, comment = "#"))
    #TO DO: fix this so it will join on biopsy_id or sample_id depending on which one is present, Done?
    sample_table = left_join(sample_table, this_data)
  }
  return(sample_table)
}
