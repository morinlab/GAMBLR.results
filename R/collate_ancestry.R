#' @title Collate Ancestry.
#'
#' @description Gather ancestry information and expand the incoming sample table (or metadata).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param somalier_output Somalier ancestery.tsv
#'
#' @return A table.
#'
#' @import stringr readr dplyr
#'
#' @noRd
#'
#' @examples
#' table = collate_ancestry(sample_table = "my_sample_table.txt")
#'
collate_ancestry = function(sample_table,
                            seq_type_filter="genome",
                            somalier_output){
  if(seq_type_filter=="capture"){
    message("skipping ancestry for this seq_type")
    return(sample_table)
  }
  if(missing(somalier_output)){
    somalier_output = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/gambl/somalier_current/02-ancestry/2020_08_07.somalier-ancestry.tsv"
  }
  somalier_all = suppressMessages(read_tsv(somalier_output))
  somalier_all = mutate(somalier_all, sample_id = str_remove(`#sample_id`, pattern = ":.+")) %>%
    dplyr::select(-`#sample_id`, -given_ancestry)
  somalier_all = dplyr::select(somalier_all, sample_id, predicted_ancestry, PC1, PC2, PC3, PC4, PC5)
  sample_table = left_join(sample_table, somalier_all)
  return(sample_table)
}
