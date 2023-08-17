#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range.
#'
#' @details This function returns CN segments for s specified region.
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
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param with_chr_prefix Boolean parameter for toggling if chr prefixes should be present in the return, default is FALSE.
#' @param streamlined Return a basic rather than full MAF format. Default is FALSE.
#' @param from_flatfile Set to TRUE by default.
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr readr RMariaDB DBI glue
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
                           chromosome,
                           qstart,
                           qend,
                           projection = "grch37",
                           this_seq_type = "genome",
                           with_chr_prefix = FALSE,
                           streamlined = FALSE,
                           from_flatfile = TRUE){

  #checks
  remote_session = check_remote_configuration(auto_connect = TRUE)

  #get wildcards from this_seq_type (lazy)
  seq_type = this_seq_type

  #perform wrangling on the region to have it in the correct format.
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }

  #deal with chr prefixes for region, based on selected genome projection.
  if(projection == "grch37"){
    if(grepl("chr", chromosome)){
      chromosome = gsub("chr", "", chromosome)
    }
  }else{
    if(!grepl("chr", chromosome)){
      chromosome = paste0("chr", chromosome)
    }
  }

  #enforce data type for qend and qstart coordiantes.
  qstart = as.numeric(qstart)
  qend = as.numeric(qend)

  if(from_flatfile){
    cnv_flatfile_template = check_config_value(config::get("results_flatfiles")$cnv_combined$icgc_dart)
    cnv_path =  glue::glue(cnv_flatfile_template)
    full_cnv_path =  paste0(check_config_value(config::get("project_base")), cnv_path)

    #check permissions to ICGC data.
    permissions = file.access(full_cnv_path, 4)
    if(permissions == -1){
      message("restricting to non-ICGC data")
      cnv_flatfile_template = check_config_value(config::get("results_flatfiles")$cnv_combined$gambl)
      cnv_path =  glue::glue(cnv_flatfile_template)
      full_cnv_path =  paste0(check_config_value(config::get("project_base")), cnv_path)
    }

    #check for missingness.
    if(!file.exists(full_cnv_path)){
      print(paste("missing: ", full_cnv_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    all_segs = suppressMessages(read_tsv(full_cnv_path)) %>%
      dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
      as.data.frame()

  }else{
    db = check_config_value(config::get("database_name"))
    table_name = check_config_value(config::get("results_tables")$copy_number)
    table_name_unmatched = check_config_value(config::get("results_tables")$copy_number_unmatched)
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

    all_segs_matched = dplyr::tbl(con, table_name) %>%
      dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
      as.data.frame() %>%
      dplyr::mutate(method = "battenberg")

    all_segs_unmatched = dplyr::tbl(con, table_name_unmatched) %>%
      dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
      as.data.frame() %>%
      dplyr::filter(! ID %in% all_segs_matched$ID)  %>%
      dplyr::mutate(method = "controlfreec")

      DBI::dbDisconnect(con)

      all_segs = rbind(all_segs_matched, all_segs_unmatched)
  }

  #mutate CN states.
  all_segs = dplyr::mutate(all_segs, CN = round(2*2^log.ratio))

  #deal with chr prefixes
  if(!with_chr_prefix){
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

  #subset to only a few columns with streamlined = TRUE.
  if(streamlined){
    all_segs = dplyr::select(all_segs, ID, CN)
  }

  #return data frame with CN segments
  return(all_segs)
}
