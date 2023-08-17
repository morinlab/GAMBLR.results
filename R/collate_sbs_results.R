#' @title Collate SBS Results.
#'
#' @description Bring in the results from mutational signature analysis.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param file_path Optional path to SBS file.
#' @param scale_vals Parameter not used?
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#'
#' @return A data frame with new columns added.
#'
#' @import dplyr tibble GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' collated = collate_sbs_results(sample_table = sample_table,
#'                                sbs_manipulation = sbs_manipulation)
#
collate_sbs_results = function(sample_table,
                               seq_type_filter = "genome",
                               file_path,
                               scale_vals = FALSE,
                               sbs_manipulation = ""){
  if(seq_type_filter!="genome"){
    message("skipping sbs for seq_type")
    return(sample_table)
  }
  if(missing(file_path)){
    base = GAMBLR.helpers::check_config_value(config::get("project_base"))

    file_path = paste0(base,"icgc_dart/sigprofiler-1.0/02-extract/genome--hg38/BL_HGBL_DLBCL_FL_COMFL_CLL_MCL_B-ALL_PBL_DLBCL-BL-like_UNSPECIFIED_SCBC_MM_all/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt")
  }
  message(paste("loading",file_path))
  signatures = read.csv(file_path, sep = "\t", header = 1, row.names = 1)
  rs = rowSums(signatures)
  cn = colnames(signatures)
  new_sig = signatures
  if(sbs_manipulation == "scale"){
    #for(col in cn){
    #  scaled_vals = signatures[,col] / rs
    #  new_sig[,col]=scaled_vals
    #}
    #sbs1 = signatures[,"SBS1"] / rs
    #sbs5 = signatures[,"SBS5"] / rs
    #sbs9 = signatures[,"SBS9"] / rs
    #sbs8 = signatures[,"SBS8"] / rs
    #sbs12 = signatures[,"SBS12"] / rs
    #sbs17b = signatures[,"SBS17b"]/rs
    #sbs18 = signatures[,"SBS18"]/rs
    #sbs84 = signatures[,"SBS84"]/rs
    #sbs85 = signatures[,"SBS85"]/rs
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
    sbs = apply(signatures, 2, function(x){x/rs}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
  }
  else if(sbs_manipulation == "log"){
    #sbs1 = log(signatures[,"SBS1"]+1)
    #sbs5 = log(signatures[,"SBS5"]+1)
    #sbs9 = log(signatures[,"SBS9"]+1)
    #sbs8 = log(signatures[,"SBS8"]+1)
    #sbs12 = log(signatures[,"SBS12"]+1)
    #sbs17b = log(signatures[,"SBS17b"]+1)
    #sbs18 = log(signatures[,"SBS18"]+1)
    #sbs84 = log(signatures[,"SBS84"]+1)
    #sbs85 = log(signatures[,"SBS85"]+1)
    sbs = apply(signatures, 2, function(x){log(x + 1)}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
  }else{
    #sbs1 = signatures[,"SBS1"]
    #sbs5 = signatures[,"SBS5"]
    #sbs9 = signatures[,"SBS9"]
    #sbs8 = signatures[,"SBS8"]
    #sbs12 = signatures[,"SBS12"]
    #sbs17b = signatures[,"SBS17b"]
    #sbs18 = signatures[,"SBS18"]
    #sbs84 = signatures[,"SBS84"]
    #sbs85 = signatures[,"SBS85"]
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
    sbs = signatures %>%
      rownames_to_column("sample_id")
  }
  sample_table = left_join(sample_table, sbs)
  return(sample_table)
}
