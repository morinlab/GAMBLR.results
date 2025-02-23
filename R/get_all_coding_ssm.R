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
#' my_meta = suppressMessages(get_gambl_metadata())
#' maf_all_seqtype = get_all_coding_ssm(my_meta)
#' table(maf_all_seqtype$maf_seq_type)
#'
get_all_coding_ssm = function(these_samples_metadata = NULL,
                              include_silent=FALSE,
                              ...){
  if(missing(these_samples_metadata)){
    stop("these_samples_metadata is required")
  }
  these_samples_metadata = dplyr::filter(these_samples_metadata,
                                         seq_type!="mrna")
  seq_types_in_metadata = unique(these_samples_metadata$seq_type)
  if("capture" %in% seq_types_in_metadata){
    capture_ssm = suppressWarnings(suppressMessages(get_coding_ssm(these_samples_metadata = 
                                 dplyr::filter(these_samples_metadata,
                                               seq_type=="capture"),
                               this_seq_type = "capture", ...))) 
  }
  if("genome" %in% seq_types_in_metadata){
    genome_ssm = suppressWarnings(suppressMessages(get_coding_ssm(these_samples_metadata = 
                                dplyr::filter(these_samples_metadata,
                                              seq_type=="genome"),
                              this_seq_type = "genome", ...)))
     
  }
  if(length(seq_types_in_metadata)>1){
    merged_ssm = bind_genomic_data(capture_ssm,genome_ssm)
    return(merged_ssm)
  }else if("capture" %in% seq_types_in_metadata){
    return(capture_ssm)
  }else if("genome" %in% seq_types_in_metadata){
    return(genome_ssm)
  }
}
