#' @title Get CN States.
#'
#' @description Get a copy number matrix for all samples based on segmented data in the database.
#' This function is deprecated and has been replaced by 
#' GAMBLR.utils::segmented_data_to_cn_matrix and GAMBLR.results::get_cn_segments
#'
#' @details This function returns CN states for the specified regions using the CN data
#' from GAMBLR.results and (optionally) assumes regions with no data are diploid.
#' For how to determine/specify the coordinates of each region, refer to the
#' parameter descriptions and examples.
#'
#'
#' @return Copy number matrix with sample_id as rows and regions as columns.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#'   get_cnv_and_ssm_status()
#' }
get_cn_states = function(){
  stop("this function is deprecated. Please use get_cnv_and_ssm_status or GAMBLR.utils::segmented_data_to_cn_matrix")
}
