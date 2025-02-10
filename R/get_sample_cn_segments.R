#' @title Get Sample CN Segments.
#'
#' @description Get all segments for a single (or multiple) sample_id(s).
#'
#' @details This function returns CN segments for samples. This works for single sample or multiple samples.
#' Specify the sample IDs you are interested in with `these_sample_ids` (as a vector of characters),
#' Or call this function with `these_samples_metadata` if you already have a metadata table subset to the sample IDs of interest.
#' If none off the above parameters are specified, the function will return CN segments for available samples.
#' Note, this. function internally calls helper function id_ease for dealing with sample IDs and metadata tables.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::assign_cn_to_ssm], [GAMBLR.results::get_cn_segments], [GAMBLR.results::get_cn_states],
#'
#' @param these_sample_ids Optional argument, sample_id (vector of characters) for the sample(s) to retrieve segments for. If not provided, the function will return CN segments for all available sample IDs present in the current metadata.
#' @param these_samples_metadata Optional, provide a metadata (data frame) subset to the sample IDs of interest.
#' @param from_flatfile Set to TRUE by default.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param with_chr_prefix Set to TRUE to add a chr prefix to chromosome names. Default is FALSE.
#' @param streamlined Return a minimal output rather than full details. Default is FALSE.
#'
#' @return A data frame of segments for a specific or multiple sample ID(s).
#'
#' @import dplyr readr RMariaDB DBI
#' @export
#' @keywords internal

get_sample_cn_segments = function(these_sample_ids = NULL,
                                  these_samples_metadata = NULL,
                                  from_flatfile = TRUE,
                                  projection = "grch37",
                                  this_seq_type = "genome",
                                  with_chr_prefix = FALSE,
                                  streamlined = FALSE){
  stop("this function is deprecated. Please use get_cn_segments instead")
}
