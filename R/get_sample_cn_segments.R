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
#'
#' @examples
#' #return everything (default: genome, grch37)
#' all_segs = get_sample_cn_segments()
#'
#' #using sample ID parameter
#' one_sample = get_sample_cn_segments(these_sample_ids = "00-22011_tumorB",
#'                                     this_seq_type = "capture")
#'
#' two_samples = get_sample_cn_segments(these_sample_ids = c("00-14595_tumorA",
#'                                                           "00-16220_tumorB"),
#'                                      projection = "hg38")
#'
#' #get some metadata
#' this_meta = get_gambl_metadata() %>%
#'  dplyr::filter(sample_id %in% c("BLGSP-71-23-00440-01A-01E",
#'                                 "FL1014T2",
#'                                 "HTMCP-01-06-00563-01A-01D"))
#'
#' capture_meta = get_gambl_metadata(seq_type_filter = "capture")
#'
#' #return sample segments by using the metadata parameter
#' segs_meta = get_sample_cn_segments(these_samples_metadata = this_meta)
#'
#' all_segs_cap = get_sample_cn_segments(these_samples_metadata = capture_meta,
#'                                       projection = "hg38")
#'
get_sample_cn_segments = function(these_sample_ids = NULL,
                                  these_samples_metadata = NULL,
                                  from_flatfile = TRUE,
                                  projection = "grch37",
                                  this_seq_type = "genome",
                                  with_chr_prefix = FALSE,
                                  streamlined = FALSE){
  #get sample IDs
  meta_ids = id_ease(these_sample_ids = these_sample_ids,
                     these_samples_metadata = these_samples_metadata,
                     this_seq_type = this_seq_type)

  #subset returned list to sample IDs
  these_samples = meta_ids$sample_id

  if(from_flatfile){
    seq_type = this_seq_type
    cnv_flatfile_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$cnv_combined$icgc_dart)
    cnv_path =  glue::glue(cnv_flatfile_template)
    full_cnv_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), cnv_path)
    local_full_cnv_path =  paste0(config::get("project_base"), cnv_path)

    if(file.exists(local_full_cnv_path)){
      full_cnv_path = local_full_cnv_path
      #use local file when available
    }

    # check permissions to ICGC data
    permissions = file.access(full_cnv_path, 4)
    if(permissions == -1){
      message("restricting to non-ICGC data")
      cnv_flatfile_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$cnv_combined$gambl)
      cnv_path =  glue::glue(cnv_flatfile_template)
      full_cnv_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), cnv_path)
    }

    #check for missingness
    if(!file.exists(full_cnv_path)){
      print(paste("missing: ", full_cnv_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    #read seg file and subset to sample IDs of interest
    all_segs = suppressMessages(read_tsv(full_cnv_path))
    all_segs = dplyr::filter(all_segs, ID %in% these_samples)

  }else{
    stop("Database support is deprecated, please use this function with `from_flatfile = TRUE`")
  }

  all_segs = dplyr::mutate(all_segs, CN = round(2*2^log.ratio))

  #deal with chr prefixes
  if(!with_chr_prefix){
    all_segs = all_segs %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom))
  }else{
    if(!grepl("chr", all_segs$chrom[1])){
      all_segs$chrom = paste0("chr", all_segs$chrom)
    }
  }

  if(streamlined){all_segs = dplyr::select(all_segs, ID, CN)}

  return(all_segs)
}
