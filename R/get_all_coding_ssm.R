#' @title Get all Coding SSMs
#'
#' @description Retrieve all coding SSMs from GAMBL in MAF-like format,
#' regardless of seq_type.
#'
#' @details Effectively retrieve coding SSM calls from one or all DNA seq_type.
#' For additional optional arguments, see [GAMBLR.results::get_coding_ssm]

#' @param these_samples_metadata Supply a metadata table containing the
#' sample/seq_type combinations you want. 
#' @param include_silent If set to TRUE, silent/synonymous mutations in the
#' coding regions will also be returned. Default is FALSE.
#' @param ... Additional arguments passed to get_coding_ssm
#'
#' @param with_metadata Set to TRUE to return a list with the mutations
#' (`maf_data`) and `metadata`: `these_samples_metadata`, all rows kept,
#' with an `ssm_available` column (and `all_available`) saying which samples
#' have SSM data, so samples without mutations can be told from samples
#' without data. See [annotate_data_availability]. Default FALSE.
#' @return A data frame containing all the MAF data columns (one row per
#' mutation).
#'
#' @import dplyr tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' my_meta = suppressMessages(get_gambl_metadata())
#' maf_all_seqtype = get_all_coding_ssm(my_meta)
#' 
#' table(maf_all_seqtype$maf_seq_type)
#' 
#' # most common mutations by gene and Variant_Classification
#' dplyr::group_by(maf_all_seqtype,Hugo_Symbol,Variant_Classification) %>% 
#'   dplyr::count() %>% 
#'   dplyr::arrange(desc(n))
get_all_coding_ssm = function(these_samples_metadata = NULL,
                              include_silent=FALSE,
                              with_metadata = FALSE,
                              ...){
  if(with_metadata){
    if(missing(these_samples_metadata) || is.null(these_samples_metadata)){
      stop("with_metadata = TRUE needs these_samples_metadata")
    }
    return(.with_metadata(match.call(), parent.frame(), "maf_data", function()
      annotate_ssm_availability(these_samples_metadata, source = "cds",
                                projection = if(is.null(list(...)$projection)) "grch37" else list(...)$projection,
                                augmented = !isFALSE(list(...)$augmented))))
  }
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
                               this_seq_type = "capture",
                               include_silent = include_silent,
                               ...))) 
  }
  if("genome" %in% seq_types_in_metadata){
    genome_ssm = suppressWarnings(suppressMessages(get_coding_ssm(these_samples_metadata = 
                                dplyr::filter(these_samples_metadata,
                                              seq_type=="genome"),
                              this_seq_type = "genome", 
                              include_silent = include_silent,
                              ...)))
     
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
