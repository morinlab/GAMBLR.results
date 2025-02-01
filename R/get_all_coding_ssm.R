#' @title Get all Coding SSMs
#'
#' @description Retrieve all coding SSMs from GAMBL in MAF-like format, regardless of seq_type.
#'
#' @details Effectively retrieve coding SSM calls from one or all DNA seq_type. For additional optional arguments, see [GAMBLR.results::get_coding_ssm]

#' @param these_samples_metadata Supply a metadata table containing the sample/seq_type combinations you want. 
#' @param include_silent If set to TRUE, silent/synonymous mutations in the coding regions will also be returned. 
#' @param ... Additional arguments passed to get_coding_ssm
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#'
#' @import dplyr tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' 
#' maf_data_all_exome_and_genome = get_coding_ssm(these_samples_metadata = get_gambl_metadata(seq_type_filter=c("genome","capture")))
#'
#'
get_all_coding_ssm = function(these_samples_metadata = NULL,
                              include_silent=FALSE,
                              ...){
  these_samples_metadata = dplyr::filter(these_samples_metadata,seq_type!="mrna")
  capture_ssm = get_coding_ssm(these_samples_metadata = 
                                 dplyr::filter(these_samples_metadata,
                                               seq_type=="capture"),
                               this_seq_type = "capture", ...) 

  genome_ssm = get_coding_ssm(these_samples_metadata = 
                                dplyr::filter(these_samples_metadata,
                                              seq_type=="genome"),
                              this_seq_type = "genome", ...) 

  merged_ssm = bind_genomic_data(capture_ssm,genome_ssm)
  
  return(merged_ssm)
}
