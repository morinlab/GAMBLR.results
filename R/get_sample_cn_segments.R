#' @title Get Sample CN Segments.
#'
#' @description Get all segments for a single (or multiple) sample_id(s).
#'
#' @details This function returns CN segments for samples. This works for single sample or multiple samples.
#' For multiple samples, remember to set the Boolean parameter `multiple_samples = TRUE` and give the `sample_lsit` a vector of characters with one sample ID per row.
#' For more information on how this function can be run in different ways, refer to the parameter descriptions, examples and vignettes.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::assign_cn_to_ssm], [GAMBLR.results::get_cn_segments], [GAMBLR.results::get_cn_states],
#'
#' @param this_sample_id Optional argument, single sample_id for the sample to retrieve segments for.
#' @param multiple_samples Set to TRUE to return cn segments for multiple samples specified in `samples_list` parameter. Default is FALSE.
#' @param sample_list Optional vector of type character with one sample per row, required if multiple_samples is set to TRUE.
#' @param from_flatfile Set to TRUE by default.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param with_chr_prefix Set to TRUE to add a chr prefix to chromosome names. Default is FALSE.
#' @param streamlined Return a minimal output rather than full details. Default is FALSE.
#'
#' @return A data frame of segments for a specific or multiple sample ID(s).
#'
#' @import dplyr readr RMariaDB DBI glue GAMBLR.helpers
#' @export
#'
#' @examples
#' # Return cn segments for multiple samples (read csv with one sample per line):
#' sample_list = readLines("../samples-test.csv")
#' multiple_samples = get_sample_cn_segments(multiple_samples = TRUE, sample_list = sample_list)
#' #Return cn segments for multiple samples (provided as vector of sample IDs):
#' these_sample_list = c("00-15201_tumorA", "00-15201_tumorB")
#'
#' samples = get_sample_cn_segments(multiple_samples = TRUE,
#'                                  sample_list = these_sample_list)
#' # For capture
#' samples = get_sample_cn_segments(
#'  multiple_samples = TRUE,
#'  sample_list = these_sample_list,
#'  this_seq_type = "capture"
#' )
#'
get_sample_cn_segments = function(this_sample_id,
                                  multiple_samples = FALSE,
                                  sample_list,
                                  from_flatfile = TRUE,
                                  projection = "grch37",
                                  this_seq_type = "genome",
                                  with_chr_prefix = FALSE,
                                  streamlined = FALSE){
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
    if (permissions == -1) {
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

    all_segs = suppressMessages(read_tsv(full_cnv_path))
    if (!missing(this_sample_id) & !multiple_samples) {
      all_segs = dplyr::filter(all_segs, ID %in% this_sample_id)
    } else if (!missing(sample_list)) {
      all_segs = dplyr::filter(all_segs, ID %in% sample_list)
    }

  } else {

    if(!missing(this_sample_id) & !multiple_samples){
      sample_status = get_gambl_metadata() %>%
        dplyr::filter(sample_id == this_sample_id) %>%
        pull(pairing_status)

      db = GAMBLR.helpers::check_config_value(config::get("database_name"))
      table_name = GAMBLR.helpers::check_config_value(config::get("results_tables")$copy_number)
      table_name_unmatched = GAMBLR.helpers::check_config_value(config::get("results_tables")$copy_number_unmatched)
      con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

      all_segs_matched = dplyr::tbl(con, table_name) %>%
        dplyr::filter(ID == this_sample_id) %>%
        as.data.frame() %>%
        dplyr::mutate(method = "battenberg")

      all_segs_unmatched = dplyr::tbl(con, table_name_unmatched) %>%
        dplyr::filter(ID == this_sample_id) %>%
        as.data.frame() %>%
        dplyr::filter(! ID %in% all_segs_matched$ID) %>%
        dplyr::mutate(method = "controlfreec")
    }

    if(multiple_samples & missing(this_sample_id)){
      sample_status = get_gambl_metadata() %>%
        dplyr::filter(sample_id %in% sample_list) %>%
        pull(pairing_status)

      db = GAMBLR.helpers::check_config_value(config::get("database_name"))
      table_name = GAMBLR.helpers::check_config_value(config::get("results_tables")$copy_number)
      table_name_unmatched = GAMBLR.helpers::check_config_value(config::get("results_tables")$copy_number_unmatched)
      con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

      all_segs_matched = dplyr::tbl(con, table_name) %>%
        dplyr::filter(ID %in% sample_list) %>%
        as.data.frame() %>%
        dplyr::mutate(method = "battenberg")

      all_segs_unmatched = dplyr::tbl(con, table_name_unmatched) %>%
        dplyr::filter(ID %in% sample_list) %>%
        as.data.frame() %>%
        dplyr::filter(! ID %in% all_segs_matched$ID) %>%
        dplyr::mutate(method = "controlfreec")
    }

    all_segs = rbind(all_segs_matched, all_segs_unmatched)

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
