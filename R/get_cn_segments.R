


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
#' @param max_CN Cap CN values at this value. Default is 20. This is to avoid extreme CN values dominating the colour scale in visualisations.
#' Available options are "combined" (for the merge of the best tool for each data type),
#' or one of "purecn_cnvkit", "purecn_denovo", or "battenberg".
#' Other thabn "combined", the other flavours are incomplete (e.g. limited to genome,
#' matched or capture samples).
#' @param projection Desired genome coordinate system for returned CN segments. Default is "grch37".
#' @param fill_missing_with Deprecated.
#' This argument is no longer available as we are updating how GAMBL reports CN and log.ratio values.
#' For simplicity, dummy segments (if detected) are now dropped.
#' @param adjust_for_ploidy Deprecated.
#' This argument is no longer available as we are updating how GAMBL reports CN and log.ratio values.
#'
#' @param verbose Set to TRUE for a chattier experience
#' @param this_seq_type Deprecated.
#'
#' @return A data frame with CN segments for the specified region. Includes
#' a `dummy_segment` column (0/1) flagging real segments that land on
#' coordinates commonly used as a dummy/filler region elsewhere in the
#' cohort (shared by more than one sample); the dummy/filler rows
#' themselves are dropped.
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
#' # Note the default projection is "grch37"
#'
#' genome_metadata <- suppressMessages(get_gambl_metadata()) %>%
#'   dplyr::filter(seq_type == "genome")
#'
#' # Create a metadata table with a mix of seq_types
#' mixed_seq_type_meta <- dplyr::bind_rows(capture_metadata, genome_metadata)
#'
#' ## We can load the copy number segments for all samples across seq_types
#' ## To differentiate samples that share a sample_id but differ in seq_type,
#' ## a new column seg_seq_type is added to the output
#'
#' all_seq_type_segs <- get_cn_segments(
#'   these_samples_metadata = mixed_seq_type_meta
#' )
#'
#' dplyr::group_by(all_seq_type_segs,seg_seq_type) %>% dplyr::count()
#'

get_cn_segments <- function(these_samples_metadata,
                            projection = "grch37",
                            flavour = "combined",
                            this_seq_type,
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
    if(verbose){
      print(seq_type)
      print(head(seg))
    }
    
    df_list[[seq_type]] <- seg
  }
  if (any(unique(df_list[["capture"]]$ID) %in% unique(df_list[["genome"]]$ID))) {
    stop("overlapping IDs found!")
  }
  if(verbose){
    print(names(df_list))
  }

  all_segs <- do.call("bind_rows", df_list)
  if(!"log.ratio" %in% colnames(all_segs)){
    all_segs = rename(all_segs,
      c("log.ratio"="seg.mean"))
  }
  if(!"CN" %in% colnames(all_segs) || any(is.na(all_segs[["CN"]]))){
    all_segs = dplyr::mutate(all_segs, CN = 2 * 2^log.ratio)
  }
  if(max_CN > 0){
    all_segs = dplyr::mutate(all_segs, CN = ifelse(CN > max_CN, max_CN, CN))
  }
  # If dummy segments are present, drop the placeholder rows themselves, but
  # keep a dummy_segment flag on the real rows that land on coordinates
  # commonly used as a dummy/filler region elsewhere in the cohort (shared
  # by >1 sample) -- a QC signal for low-confidence/boundary-artifact
  # regions, re-attached by coordinate rather than carried through on the
  # original (now-dropped) rows. This change could impact the GISTIC
  # functionality in GAMBLR so this may need to be revisited in the future.
  # If you want the raw dummy rows instead, load the raw seg files yourself
  # and specify fill_missing_with or adjust_for_ploidy as needed.
  if("dummy_segment" %in% colnames(all_segs)){
    if(verbose){
      message("dropping all rows where dummy_segment == 1")
    }
    dummy_segs = dplyr::filter(all_segs, dummy_segment==1)
    if(verbose){
      message(sprintf(
        "found %d dummy segment row(s) across %d sample(s)",
        nrow(dummy_segs), n_distinct(dummy_segs$ID)
      ))
    }
    # identify common dummy segment coordinates across samples and seq_types
    common_dummy_coords = dplyr::group_by(dummy_segs, chrom, start, end) %>%
      dplyr::summarise(n_samples = n_distinct(ID), .groups = "drop") %>%
      dplyr::filter(n_samples > 1) %>%
      dplyr::mutate(match_type = "exact")
    if(verbose){
      message(sprintf(
        "%d distinct dummy coordinate(s) shared by >1 sample (candidates for flagging)",
        nrow(common_dummy_coords)
      ))
    }
    # tolerate two boundary conventions a real segment's coordinates can
    # differ from the dummy coordinate by: the end alone off by 1bp
    # (end_shifted), or both the start and the end off by 1bp together
    # (start_shifted) -- confirmed against real data to be the more common
    # of the two non-exact cases by a wide margin (hundreds of rows vs.
    # effectively none for a start-only shift), so this isn't a symmetric
    # "either boundary independently" tolerance.
    right_shifted_dummy_coords = dplyr::mutate(common_dummy_coords,
      end = end + 1, match_type = "end_shifted")
    left_shifted_dummy_coords = dplyr::mutate(common_dummy_coords,
      start = start - 1, end = end + 1, match_type = "start_shifted")

    common_dummy_coords = bind_rows(common_dummy_coords, right_shifted_dummy_coords, left_shifted_dummy_coords) %>%
      # dedupe on coordinate alone, not the whole row: two different common-
      # dummy regions that are only ~1bp apart can produce the same shifted
      # (chrom, start, end) with different n_samples: left as-is, that
      # leaves duplicate join keys below, which fans out (duplicates) any
      # real segment landing on the colliding coordinate. Binding "exact"
      # first means a coordinate reachable both ways keeps the "exact"
      # label (distinct() keeps the first occurrence of each key).
      dplyr::distinct(chrom, start, end, .keep_all = TRUE)

    all_segs = dplyr::filter(all_segs, dummy_segment!=1 | is.na(dummy_segment)) %>%
      dplyr::select(-dummy_segment)
    # annotate rows that match coordinates of known dummy segments exactly
    if(nrow(common_dummy_coords) > 0){
      all_segs = dplyr::left_join(all_segs, common_dummy_coords, by = c("chrom", "start", "end")) %>%
        dplyr::mutate(dummy_segment = ifelse(!is.na(n_samples), 1, 0))
      if(verbose){
        flagged_counts = dplyr::filter(all_segs, dummy_segment == 1) %>%
          dplyr::count(match_type, name = "n_rows")
        total_flagged = sum(flagged_counts$n_rows)
        message(sprintf(
          "flagged %d real segment row(s) as dummy_segment == 1, by matching coordinate set:",
          total_flagged
        ))
        for(mt in c("exact", "end_shifted", "start_shifted")){
          n = flagged_counts$n_rows[flagged_counts$match_type == mt]
          if(length(n) == 0) n = 0
          message(sprintf("  %s: %d row(s)", mt, n))
        }
      }
      all_segs = dplyr::select(all_segs, -n_samples, -match_type)
    } else {
      all_segs = dplyr::mutate(all_segs, dummy_segment = 0)
      if(verbose){
        message("no common dummy coordinates found; dummy_segment set to 0 for all rows")
      }
    }
  }
  # return S3 class with CN segments and genome_build
  all_segs <- create_seg_data(all_segs, projection)
  return(all_segs)
}
