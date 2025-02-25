#' @title Get Sample CN Segments.
#'
#' @description Get all segments for a single (or multiple) sample_id(s).
#'
#' @details Deprecated. See [GAMBLR.results::get_cn_segments].
#' 
#' @return A data frame of segments for a specific or multiple sample ID(s).
#'
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
