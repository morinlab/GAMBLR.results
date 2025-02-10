#' Create Segmented Data Object
#'
#' This function constructs an S3 object for segmented genomic data by augmenting a
#' standard data frame with a custom class and a genome build attribute. The returned
#' object inherits from \code{data.frame} and is assigned the class \code{"seg_data"}.
#'
#' @param seg_df A data frame containing segmented data. This data frame should include
#'   the necessary columns that describe segments (for example, chromosome, start, end, and
#'   segment mean values). There is no formal check on column names; however, downstream
#'   methods may expect certain column names.
#' @param genome_build A character string specifying the genome build. Valid values are
#'   \code{"grch37"} and \code{"hg38"}. An error is thrown if an invalid genome build is provided.
#'
#' @return An object of class \code{"seg_data"} (in addition to its inherited classes) with
#'   an attached attribute \code{genome_build} containing the specified genome build.
#'
#' @details The \code{create_seg_data} function is designed to encapsulate segmented data along
#'   with its metadata, making it easier to work with such objects in downstream analyses.
#'   The function performs basic checks to ensure that the input data is a data frame and that
#'   the genome build is valid. Further processing or validations (e.g., ensuring the presence
#'   of specific columns) should be handled prior to calling this constructor if needed.
#'
#' @examples
#' \dontrun{
#' # Example segmented data
#' seg_df <- data.frame(
#'   chrom = c("chr1", "chr1", "chr2"),
#'   start = c(100000, 200000, 150000),
#'   end = c(150000, 250000, 200000),
#'   seg_mean = c(0.5, -0.3, 0.2)
#' )
#'
#' # Create a seg_data object using the hg38 genome build
#' seg_obj <- create_seg_data(seg_df, genome_build = "hg38")
#'
#' # Print the object and check its genome build attribute
#' print(seg_obj)
#' attr(seg_obj, "genome_build")
#' }
#'
#' @export
#' @keywords internal
create_seg_data <- function(seg_df, genome_build) {
  if (!inherits(seg_df, "data.frame"))
    stop("data must be a data frame")
  if (!genome_build %in% c("grch37", "hg38"))
    stop("Invalid genome build")
  #ensure chr prefixes are there when necessary 
  if(genome_build=="grch37"){
    if(all(str_detect(seg_df$chrom, "chr"))){
      seg_df = seg_df %>%
        dplyr::mutate(chrom = gsub("chr", "", chrom))
    }
  }else{
    if(all(!str_detect(seg_df$chrom, "chr"))){
      seg_df = seg_df %>%
        dplyr::mutate(chrom = paste0("chr", chrom))
    }
  }
  structure(seg_df,
            class = c("seg_data", class(seg_df)),
            genome_build = genome_build)
}

#' @export
print.seg_data <- function(x, ...) {
  cat("SEG Data Object\n")
  cat("Genome Build:", attr(x, "genome_build"), "\n")
  cat("Showing first 10 rows:\n")
  # Convert to a plain data.frame (if not already) so that printing uses the default
  # data.frame print method rather than printing as a list.
  print(utils::head(as.data.frame(x), 10))
}


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
