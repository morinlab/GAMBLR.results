#' @title Collate CSR Results
#'
#' @description Collate a few CSR annotations, including MiXCR.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return The sample table with additional columns.
#'
#' @import readr dplyr
#'
#' @noRd
#'
#' @examples
#' gambl_results_derived = collate_csr_results(gambl_results_derived)
#'
collate_csr_results = function(sample_table,
                               seq_type_filter = "genome"){

   if(seq_type_filter=="capture"){
     return(sample_table) #result doesn't exist or make sense for this seq_type
   }
   csr = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-nthomas/results/icgc_dart/mixcr_current/level_3/mixcr_genome_CSR_results.tsv"))
   sm_join = inner_join(sample_table, csr, by = c("sample_id" = "sample"))
   pt_join = inner_join(sample_table, csr, by = c("patient_id" = "sample"))
   complete_join = bind_rows(pt_join, sm_join) %>%
     bind_rows(dplyr::filter(sample_table, !patient_id %in% c(pt_join$patient_id, sm_join$patient_id))) %>%
     unique()

  return(complete_join)
}
