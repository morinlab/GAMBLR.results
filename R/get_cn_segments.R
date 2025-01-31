# Constructor function for segmented data
#' @export
create_seg_data <- function(seg_df, genome_build) {
  if (!inherits(seg_df, "data.frame")) stop("data must be a data frame")
  if (!genome_build %in% c("grch37", "hg38")) stop("Invalid genome build")
  
  structure(seg_df, 
            class = c("seg_data", class(seg_df)), 
            genome_build = genome_build)
}


#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range.
#'
#' @details This function returns CN segments for s specified region, chromosome, or with no filtering on region (i.e. everything).
#' There are multiple ways a region can be specified.
#' For example, the user can provide the full region in a "region" format (chr:start-end) to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually with `chr`, `start`, and `end` parameters.
#' For more usage examples, refer to the parameter descriptions and examples in the vignettes.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::assign_cn_to_ssm], [GAMBLR.results::get_cn_states], [GAMBLR.results::get_sample_cn_segments]
#'
#' @param region Region formatted like chrX:1234-5678 or X:1234-56789.
#' @param chromosome The chromosome you are restricting to. Required parameter if region is not specified.
#' @param qstart Start coordinate of the range you are restricting to. Required parameter if region is not specified.
#' @param qend End coordinate of the range you are restricting to. Required parameter if region is not specified.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param weighted_average Enable a more accurate estimation of the copy number for a region that uses a weighted average of all overlapping/contained segments.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param with_chr_prefix Boolean parameter for toggling if chr prefixes should be present in the return, default is FALSE.
#' @param streamlined Return a basic rather than full MAF format. Default is FALSE.
#' @param from_flatfile Set to TRUE by default.
#' @param these_samples_metadata Provide a metadata table to restrict the data to the samples in your table
#' @param seg_data Optionally provide the function with a data frame of segments that will be used instead of the GAMBL flatfiles
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr readr RMariaDB DBI glue GAMBLR.helpers
#' @export
#'
#' @examples
#' #Example using chromosome, qstart and qend parameters:
#' segments_region_grch37 = get_cn_segments(chromosome = "chr8",
#'                                          qstart = 128723128,
#'                                          qend = 128774067)
#' # Example for the capture samples:
#' capture_segments_region_grch37 = get_cn_segments(
#'  chromosome = "chr8",
#'  qstart = 128723128,
#'  qend = 128774067,
#'  this_seq_type = "capture"
#' )
#'
#' #Example using the regions parameter:
#' segments_region_hg38 = get_cn_segments(region = "chr8:128,723,128-128,774,067",
#'                                        projection = "hg38",
#'                                        with_chr_prefix = TRUE)
#'
get_cn_segments = function(region,
                           projection = "grch37",
                           weighted_average=FALSE,
                           seg_data,
                           this_seq_type,
                           with_chr_prefix = FALSE,
                           streamlined = FALSE,
                           from_flatfile = TRUE,
                           these_samples_metadata){

  if(!missing(seg_data) & !missing(region)){
    #for backwards compatability, this function can still do the processing of pre-loaded seg_data
    print("calling helper function")
    seg_data = process_cn_segments_by_region(seg_data=seg_data,
                                             region=region,
                                             streamlined=streamlined,
                                             weighted_average=weighted_average)
    return(seg_data)
    
  }
  if(!missing(this_seq_type)){
    message("this_seq_type has been deprecated in get_cn_segments")
  }
  if(missing(these_samples_metadata)){
      message("no metadata provided, will get segments for every genome and capture sample")
      these_samples_metadata = get_gambl_metadata() %>% filter(seq_type %in% c("genome","capture"))
      seq_type = pull(these_samples_metadata,seq_type) %>% unique()
  }else if(!missing(these_samples_metadata)){
    seq_type = pull(these_samples_metadata,seq_type) %>% unique()
  }
  cnv_flatfile_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$cnv_combined$icgc_dart)
  cnv_path =  glue::glue(cnv_flatfile_template)
  full_cnv_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), cnv_path)
  print(full_cnv_path)
    #check permissions to ICGC data.
  permissions = file.access(full_cnv_path[1], 4)
  if(permissions == -1){
    message(paste("failed loading from",full_cnv_path[1]))
    message("restricting to non-ICGC data")
    cnv_flatfile_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$cnv_combined$gambl)
    cnv_path =  glue::glue(cnv_flatfile_template)
    full_cnv_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), cnv_path)
  }

  #check for missingness.
  if(!file.exists(full_cnv_path[1])){
    print(paste("missing: ", full_cnv_path[1]))
    message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
    message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
  }
  df_list <- suppressMessages(map(full_cnv_path, read_tsv))
  
  all_segs = bind_rows(df_list) %>%  as.data.frame()
  
  all_segs = dplyr::mutate(all_segs, CN = round(2*2^log.ratio))
  
  #ensure chr prefixes are there when necessary 
  if(projection=="grch37"){
    if(all(str_detect(all_segs$chrom, "chr"))){
      all_segs = all_segs %>%
        dplyr::mutate(chrom = gsub("chr", "", chrom))
    }
  }else{
    if(all(!str_detect(all_segs$chrom, "chr"))){
      all_segs = all_segs %>%
        dplyr::mutate(chrom = paste0("chr", chrom))
    }
  }


  if(!missing(these_samples_metadata)){
    these_samples = pull(these_samples_metadata,sample_id)
    all_segs = dplyr::filter(all_segs,ID %in% these_samples)
  }

  #return S3 class with CN segments and genome_build 
  all_segs = create_seg_data(all_segs,projection)
  return(all_segs)
}
