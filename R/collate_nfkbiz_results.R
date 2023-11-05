#' @title Collate NFKBIZ Results.
#'
#' @description Determine which cases have NFKBIZ UTR mutations.
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
collate_nfkbiz_results = function(sample_table,
                                  seq_type_filter = "genome"){

  #TO DO: Update to work with hg38 projection
  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter=seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  this_region = "chr3:101578214-101578365"
  nfkbiz_ssm = get_ssm_by_region(region = this_region,this_seq_type = seq_type_filter) %>%
    pull(Tumor_Sample_Barcode) %>%
    unique
  if(seq_type_filter=="genome"){
    nfkbiz_sv = get_manta_sv(region = this_region) %>%
      pull(tumour_sample_id) %>%
      unique
    nfkbiz = unique(c(nfkbiz_ssm, nfkbiz_sv))
  }
  else{
    nfkbiz = unique(nfkbiz_ssm)
  }


  sample_table$NFKBIZ_UTR = "NEG"
  sample_table[sample_table$sample_id %in% nfkbiz, "NFKBIZ_UTR"] = "POS"
  return(sample_table)
}
