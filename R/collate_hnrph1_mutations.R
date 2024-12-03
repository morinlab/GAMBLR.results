#' @title Collate HNRNPH1 mutations.
#'
#' @description Determine which cases have HNRNPH1 regulatory mutations.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return Samples table.
#'
#' @import dplyr
#'
#' @noRd
#'
#' @examples
#' sample_table = collate_nfkbiz_results(sample_table = sample_table)
#'
collate_hnrnph1_mutations = function(sample_table,
                                  seq_type_filter = "genome"){
  
  #TO DO: Update to work with hg38 projection
  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter=seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  this_region = "chr5:179,046,194-179,046,508"
  hnrnph1 = get_ssm_by_region(region = this_region,this_seq_type = seq_type_filter) %>%
    pull(Tumor_Sample_Barcode) %>%
    unique
  
  
  
  sample_table$HNRNPH1_splice = "NEG"
  sample_table[sample_table$sample_id %in% hnrnph1, "HNRNPH1_splice"] = "POS"
  return(sample_table)
}
