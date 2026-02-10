


#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL outputs
#'
#' @details This merely loads and returns all the seg_data available for a projection
#' (genome build) and can assign a single value to dummy segments if they are
#' present/identified in the source file
#' @param these_samples_metadata User must provide a metadata table to restrict the data to the samples in your table.
#' The metadata also ensures the proper handling of duplicate sample_id across seq_types and ensures the
#' seq_type in the metadata faithfully represents the seq_type of the data
#' @param flavour Specify what pipeline or source of data to use.
#' Available options are "combined" (for the merge of the best tool for each data type),
#' or one of "purecn_cnvkit", "purecn_denovo", or "battenberg".
#' Other thabn "combined", the other flavours are incomplete (e.g. limited to genome,
#' matched or capture samples).
#' @param projection Desired genome coordinate system for returned CN segments. Default is "grch37".
#' @param fill_missing_with Specify how to fill values in dummy segments that were created to satisfy GISTIC.
#' The default is "nothing", which causes these to be dropped so empty regions
#' can be handled in subsequent processing steps. For creating a GISTIC input,
#' you would typically want to set this to "avg_ploidy".
#' This is taken care of for you by [GAMBLR.utils::prepare_gistic_inputs]
#' @param adjust_for_ploidy Set to TRUE to force segments to be adjusted
#' for the average genome-wide copy number (ploidy) of the sample. Beware! There are multiple
#' GAMBLR.* functions that can do this adjustment.
#' If you apply it at this stage you should *NOT* apply it again to the resulting data!
#'
#' @param verbose Set to TRUE for a chattier experience
#' @param this_seq_type Deprecated.
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr readr glue GAMBLR.utils
#' @export
#'
#' @examples
#' # Example for just exome/capture samples:
#' # Get metadata for just a few capture samples
#' capture_metadata <- suppressMessages(get_gambl_metadata()) %>%
#'   dplyr::filter(seq_type == "capture") %>%
#'   head()
#'
#' # Load the copy number segments for capture samples using hg38 projection
#' capture_segments_hg38 <- get_cn_segments(
#'   these_samples_metadata = capture_metadata,
#'   projection = "hg38"
#' )
#' print(capture_segments_hg38)
#'
#' genome_metadata <- suppressMessages(get_gambl_metadata()) %>%
#'   dplyr::filter(seq_type == "genome")
#' # Create a metadata table with a mix of seq_types
#' mixed_seq_type_meta <- dplyr::bind_rows(capture_metadata, genome_metadata)
#' ## We can load the copy number segments for all samples across seq_types
#' all_seq_type_segs <- get_cn_segments(
#'   these_samples_metadata = mixed_seq_type_meta
#' )
#'
#' dplyr::group_by(all_seq_type_segs,seg_seq_type) %>% count()
#'
#' # Note the default projection is "grch37"
get_cn_segments <- function(these_samples_metadata,
                            projection = "grch37",
                            flavour = "combined",
                            this_seq_type,
                            fill_missing_with = "nothing",
                            adjust_for_ploidy = FALSE,
                            max_CN = 20,
                            verbose = FALSE) {
  if (!missing(this_seq_type)) {
    stop("this_seq_type is deprecated. Subset your metadata instead.")
  }
  if (missing(these_samples_metadata)) {
    message("no metadata provided")
    message("will get segments for all available genome and capture samples")
    these_samples_metadata <- suppressMessages(get_gambl_metadata()) %>%
      dplyr::filter(seq_type %in% c("genome", "capture"))
    seq_types <- pull(these_samples_metadata, seq_type) %>% unique()
  } else {
    these_samples_metadata <- dplyr::filter(these_samples_metadata,
      seq_type %in% c("genome", "capture"))
    seq_types <- pull(these_samples_metadata, seq_type) %>% unique()
  }

  genome_ids <- dplyr::filter(these_samples_metadata,
    seq_type == "genome") %>%
    pull(sample_id)
  capture_ids <- dplyr::filter(these_samples_metadata,
    seq_type == "capture") %>%
    pull(sample_id)
  if (flavour == "combined") {
    cnv_flatfile_template <- check_config_and_value(
      "results_flatfiles$cnv_combined$icgc_dart"
    )
    coltypes = "cciiid"
  } else if (flavour == "battenberg") {
    seq_types = "genome"
    cnv_flatfile_template <- check_config_and_value(
      "results_merged$battenberg"
    )
    coltypes = "cciiidi"
  } else if(grepl("purecn",flavour)){
    seq_types = "capture"
    coltypes = "cciidd"
    cnv_flatfile_template <- check_config_and_value(
      paste0("results_merged$",flavour))
  } else{
    stop("currently available flavours: combined, purecn_denovo, purecn_cnvkit or battenberg")
  }
  df_list <- list()
  for (seq_type in seq_types) {
    cnv_path <- glue::glue(cnv_flatfile_template)
    full_cnv_path <- paste0(check_config_and_value("project_base"), cnv_path)
    # check permissions to ICGC data.
    permissions <- file.access(full_cnv_path, 4)
    if (permissions == -1) {
      message(paste("failed loading from", full_cnv_path[1]))
      message("restricting to non-ICGC data")
      cnv_flatfile_template <- check_config_and_value("results_flatfiles$cnv_combined$gambl")
      cnv_path <- glue::glue(cnv_flatfile_template)
      full_cnv_path <- paste0(check_config_and_value("project_base"), cnv_path)
    }

    # check for missingness.
    if (!file.exists(full_cnv_path)) {
      print(paste("missing: ", full_cnv_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    if (seq_type == "capture") {
      if(verbose){
        print(full_cnv_path)
      }

      seg <- suppressMessages(read_tsv(full_cnv_path,
                                      col_types = coltypes,
                                      na = c("NA", "NaN"),
                                      progress = FALSE)) %>%
        dplyr::filter(ID %in% capture_ids)
    } else {
      if(verbose){
        print(full_cnv_path)
      }

      seg <- suppressMessages(read_tsv(full_cnv_path,
                                       col_types = coltypes,
                                       na = c("NA", "NaN"),
                                       progress = FALSE)) %>%
        dplyr::filter(ID %in% genome_ids)
    }
    seg <- mutate(seg, seg_seq_type = seq_type)
    df_list[[seq_type]] <- seg
  }
  if (any(unique(df_list[["capture"]]$ID) %in% unique(df_list[["genome"]]$ID))) {
    stop("overlapping IDs found!")
  }

  all_segs <- do.call("bind_rows", df_list)
  if(!"log.ratio" %in% colnames(all_segs)){
    all_segs = rename(all_segs,
      c("log.ratio"="seg.mean"))
  }
  all_segs <- dplyr::mutate(all_segs, CN = ifelse(log.ratio == -2,0,2 * 2^log.ratio))



  if (adjust_for_ploidy || "dummy_segment" %in% colnames(all_segs)) {
    if(!"dummy_segment" %in% colnames(all_segs)){
      all_segs$dummy_segment = 0
    }
    if (fill_missing_with == "diploid") {
      if(verbose){
        print("Using diploid")
      }

      all_segs <- mutate(all_segs,
        CN = ifelse(dummy_segment == 1, 2, CN),
        log.ratio = ifelse(dummy_segment == 1, 0, log.ratio)
      )
    }
    if (adjust_for_ploidy || fill_missing_with == "avg_ploidy") {
      if(verbose){
        print("Calculating sample-wide Avg_ploidy")
      }

      real_segs <- dplyr::filter(
        all_segs,
        dummy_segment == 0
      )
      real_segs <- real_segs %>%
        mutate(
          length = end - start + 1,
          CN_seg = CN * length,
          logr_seg = log.ratio * length
        ) %>%
        group_by(ID) %>%
        summarise(
          mean = mean(CN),
          real_mean = sum(CN_seg) / sum(length),
          real_mean_logr = sum(logr_seg) / sum(length)
        ) # actual average per base
      all_segs <- left_join(all_segs, select(real_segs, ID, real_mean, real_mean_logr), by = "ID")
      all_segs = mutate(all_segs,CN=ifelse(CN>max_CN,max_CN,CN))
      all_segs <- mutate(all_segs,
        CN = ifelse(dummy_segment == 1, real_mean, CN),
        log.ratio = ifelse(dummy_segment == 1, real_mean_logr, log.ratio),
        adjusted_CN = 2 + CN - real_mean,
        rounded_CN = round(adjusted_CN, 0)
      )
    }
    #if (fill_missing_with == "nothing") {
    #  #drop dummy segments entirely
    #  all_segs <- dplyr::filter(all_segs, dummy_segment == 0)
    #} else {
    #  stop("fill_missing_with must be 'nothing', 'diploid', or 'avg_ploidy'")
    #}
  }
  #else{
  #  message("dummy segments are not annotated in the inputs")
  #  message("fill_missing_with parameter will be ignored")
  #}
  # return S3 class with CN segments and genome_build
  all_segs <- create_seg_data(all_segs, projection)
  return(all_segs)
}
