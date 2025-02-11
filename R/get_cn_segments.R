
#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL outputs
#'
#' @details This function merely loads and returns all the seg_data available for a projection (genome build)
#' @param these_samples_metadata User must provide a metadata table to restrict the data to the samples in your table. 
#' The metadata also ensures the proper handling of duplicate sample_id across seq_types and ensures the 
#' seq_type in the metadata faithfully represents the seq_type of the data
#' @param projection Desired genome coordinate system for returned CN segments. Default is "grch37".
#' @param this_seq_type Deprecated.
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr readr glue
#' @export
#'
#' @examples
#' # Example for the capture samples:
#' # Get metadata for just a few capture samples
#' capture_metadata = suppressMessages(get_gambl_metadata()) %>%
#'                       dplyr::filter(seq_type=="capture") %>% head()
#'
#' # Load the copy number segments for the capture samples using hg38 projection
#' capture_segments_hg38 = get_cn_segments(
#'                              these_samples_metadata = capture_metadata,
#'                              projection="hg38")
#' print(capture_segments_hg38)
#'
#' genome_metadata = suppressMessages(get_gambl_metadata()) %>%
#'                       dplyr::filter(seq_type=="genome") %>% head()
#' # Create a metadata table with a mix of seq_types
#' mixed_seq_type_meta = dplyr::bind_rows(capture_metadata,genome_metadata)
# # We can load the copy number segments for all samples across seq_types
#' capture_segments_default = get_cn_segments(
#'                              these_samples_metadata = mixed_seq_type_meta)
#' dplyr::group_by(capture_segments_default, ID) %>% 
#' dplyr::summarize(n=dplyr::n())
#' # Note the default projection is "grch37"
#' print(capture_segments_default)
get_cn_segments = function(these_samples_metadata,
                           projection = "grch37",
                           this_seq_type){

  if(!missing(this_seq_type)){
    stop("this_seq_type has been deprecated in get_cn_segments. Subset your metadata instead.")
  }
  if(missing(these_samples_metadata)){
      message("no metadata provided, will get segments for every genome and capture sample")
      these_samples_metadata = get_gambl_metadata() %>% filter(seq_type %in% c("genome","capture"))
      seq_types = pull(these_samples_metadata,seq_type) %>% unique()
  }else{
    these_samples_metadata = dplyr::filter(these_samples_metadata,seq_type %in% c("genome","capture"))
    seq_types = pull(these_samples_metadata,seq_type) %>% unique()
  }
  
  genome_ids = dplyr::filter(these_samples_metadata,seq_type=="genome") %>% pull(sample_id)
  capture_ids = dplyr::filter(these_samples_metadata,seq_type=="capture") %>% pull(sample_id)
  
  cnv_flatfile_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$cnv_combined$icgc_dart)
  df_list = list()
  for(seq_type in seq_types){
    cnv_path =  glue::glue(cnv_flatfile_template)
    full_cnv_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), cnv_path)
    
    #check permissions to ICGC data.
    permissions = file.access(full_cnv_path, 4)
    if(permissions == -1){
      message(paste("failed loading from",full_cnv_path[1]))
      message("restricting to non-ICGC data")
      cnv_flatfile_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$cnv_combined$gambl)
      cnv_path =  glue::glue(cnv_flatfile_template)
      full_cnv_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), cnv_path)
    }
    
    #check for missingness.
    if(!file.exists(full_cnv_path)){
      print(paste("missing: ", full_cnv_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    if(seq_type=="capture"){
      seg = suppressMessages(read_tsv(full_cnv_path)) %>% dplyr::filter(ID %in% capture_ids)
    }else{
      seg = suppressMessages(read_tsv(full_cnv_path)) %>% dplyr::filter(ID %in% genome_ids)
    }
    df_list[[seq_type]]=seg

  }
  if(any(unique(df_list[["capture"]]$ID) %in% unique(df_list[["genome"]]$ID))){
    stop("overlapping IDs found!")
  }

  all_segs = do.call("bind_rows",df_list)
  
  all_segs = dplyr::mutate(all_segs, CN = round(2*2^log.ratio))
  


  #return S3 class with CN segments and genome_build 
  all_segs = create_seg_data(all_segs,projection)
  return(all_segs)
}
