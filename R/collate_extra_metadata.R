#' @title Collate Extra Metadata.
#'
#' @description Gather additional metadata information and expand the incoming sample table (or metadata).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param file_path Path to extra metadata.
#'
#' @return A table.
#'
#' @import readr dplyr
#'
#' @noRd
#'
#' @examples
#' table = collate_extra_metadata(sample_table = "my_sample_table.txt")
#'
collate_extra_metadata = function(sample_table,
                                  file_path){

  file_path = "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/experiments/2021-04-21-Trios-MiXCR/trios_relapse_timing.tsv"
  extra_df = suppressMessages(read_tsv(file_path))
  sample_table = left_join(sample_table, extra_df, by = c("sample_id" = "biopsy_id"))
}
