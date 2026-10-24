


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
#' @param region Optional region ("chrom:start-end"): only segments
#' overlapping it are read in (saving memory). Note that whole-genome
#' segments are needed for anything that estimates a sample's baseline
#' (modal) copy number, such as GAMBLR.viz::plot_cohort_CN_region; with
#' `region`, the segments get a genome-wide `baseline_cn` (see
#' `baseline_method`).
#' @param baseline_method Add each sample's genome-wide baseline copy number
#' as `baseline_cn` (and `baseline_method`): "mode" (default; continuous),
#' "mode_rounded" (a whole number) or "mean" (see
#' [GAMBLR.utils::add_segment_baseline]); NULL adds none. These are what
#' e.g. [GAMBLR.utils::center_segments] and
#' GAMBLR.viz::plot_cohort_CN_region use. Baselines come from a cache built
#' once per seg file (see [get_cn_baselines]; the first call takes a little
#' while). Samples that can't get a reliable one (segments covering less
#' than half of the autosomal genome or fewer than 20 autosomes,
#' coordinates past the chromosome ends, or no CN) get NA, with a message;
#' functions that need a baseline stop on them, so remove them first, e.g.
#' with [annotate_cn_availability] and `cn_available`. Only for flavour
#' "combined" by default.
#' @param overlap_action What to do if a sample (profile, i.e. sample_id and
#' seq_type) has overlapping segments, which should never happen and would
#' cause double counting downstream: "warn" (default, naming the samples),
#' "error", "resolve" (clean up with
#' [GAMBLR.utils::resolve_segment_overlaps]) or "none".
#' @param this_seq_type Deprecated.
#'
#' @param with_metadata Set to TRUE to return a list with the segments
#' (`seg_data`) and `metadata`: `these_samples_metadata`, all rows kept,
#' with a `cn_available` column (and `all_available`) saying which samples
#' have usable copy number data. See [annotate_data_availability]. Default
#' FALSE.
#' @return A data frame with CN segments for the specified region. Dummy
#' (placeholder) segments are dropped, and the `dummy_segment` column with
#' them. Includes `module` (the pipeline each segment came from, e.g.
#' battenberg, cnvkit, sequenza) when the seg file has it. `CN` and
#' `log.ratio` are as in the file (CN is NA where the file has none);
#' `CN_inferred` is the copy number implied by the log ratio alone
#' (`2 * 2^log.ratio`, capped at `max_CN`), which ignores purity and ploidy
#' and so is not equivalent to `CN`.
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
                            verbose = FALSE,
                            overlap_action = "warn",
                            region = NULL,
                            baseline_method = "mode",
                            with_metadata = FALSE) {
  if(with_metadata){
    if(missing(these_samples_metadata) || is.null(these_samples_metadata)){
      stop("with_metadata = TRUE needs these_samples_metadata")
    }
    if(flavour != "combined") stop("with_metadata = TRUE is available for flavour = \"combined\"")
    return(.with_metadata(match.call(), parent.frame(), "seg_data", function()
      annotate_cn_availability(these_samples_metadata, projection = projection)))
  }
  if (!missing(this_seq_type)) {
    stop("this_seq_type is deprecated. Subset your metadata instead.")
  }
  # optional baselines: each sample's genome-wide baseline copy number, from
  # the cached per-file table (see get_cn_baselines)
  add_baselines = function(x, path) x
  # baselines come from the combined seg files' cache; other flavours get
  # them only when asked for
  if (missing(baseline_method) && flavour != "combined") baseline_method <- NULL
  if (!is.null(baseline_method)) {
    baseline_method <- match.arg(baseline_method, GAMBLR.utils::baseline_methods())
    add_baselines = function(x, path) {
      b <- .cn_baselines_for_file(path, projection = projection, verbose = verbose)
      # samples that can't get a reliable baseline (as in get_cn_baselines:
      # too sparse, coordinates past the chromosome ends, or no CN) get NA
      # rather than stopping the call; functions that need a baseline stop
      # on them, naming them
      b <- b[b$ID %in% unique(x$ID), , drop = FALSE]
      bad <- (!b$fixed_baseline %in% TRUE & (b$coverage < 0.5 | b$n_autosomes < 20)) |
        b$n_beyond_chrom_end > 0 | is.na(b$baseline_mode)
      b[bad, paste0("baseline_", GAMBLR.utils::baseline_methods())] <- NA
      if (any(bad)) {
        listed <- if (sum(bad) <= 30) b$ID[bad] else c(utils::head(b$ID[bad], 30), "...")
        message(sum(bad), " sample(s) have no baseline_cn (segments too sparse, past the ",
                "chromosome ends, or without CN): ", paste(listed, collapse = ", "),
                ". Functions that need a baseline (e.g. GAMBLR.utils::center_segments) will ",
                "stop on them; remove them with annotate_cn_availability() and cn_available.")
      }
      missing_ids <- setdiff(unique(x$ID), b$ID)
      if (length(missing_ids) > 0) {
        message(length(missing_ids), " sample(s) are missing from the baseline cache and have ",
                "no baseline_cn")
      }
      x$baseline_cn <- b[[paste0("baseline_", baseline_method)]][match(x$ID, b$ID)]
      x$baseline_method <- baseline_method
      x
    }
  }
  # optional region: only segments overlapping it are read in
  in_region = function(x) x
  if (!is.null(region)) {
    clean <- gsub(",", "", region)
    parts <- regmatches(clean, regexec("^(chr)?([0-9XY]+):([0-9.e+]+)-([0-9.e+]+)$", clean))[[1]]
    if (length(parts) == 0) stop("region must look like \"chr3:187425881-187500524\"")
    r_chrom <- parts[3]; r_start <- as.numeric(parts[4]); r_end <- as.numeric(parts[5])
    in_region = function(x) x[gsub("^chr", "", x$chrom) == r_chrom & x$end > r_start &
                                x$start < r_end, , drop = FALSE]
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
    # set per file from its header (see .combined_seg_col_types)
    coltypes = NULL
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

    if (flavour == "combined") {
      coltypes <- .combined_seg_col_types(full_cnv_path)
    }
    if (seq_type == "capture") {
      if(verbose){
        print(full_cnv_path)
      }

      # read lazily and keep only these samples' segments straight away
      seg <- suppressMessages(read_tsv(full_cnv_path,
                                       col_types = coltypes,
                                       na = c("NA", "NaN"),
                                       progress = FALSE, lazy = TRUE))
      seg <- in_region(seg[seg$ID %in% capture_ids, , drop = FALSE])
      seg <- add_baselines(seg, full_cnv_path)
    } else {
      if(verbose){
        print(full_cnv_path)
      }

      seg <- suppressMessages(read_tsv(full_cnv_path,
                                       col_types = coltypes,
                                       na = c("NA", "NaN"),
                                       progress = FALSE, lazy = TRUE))
      seg <- in_region(seg[seg$ID %in% genome_ids, , drop = FALSE])
      seg <- add_baselines(seg, full_cnv_path)
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
  # CN and log.ratio are kept as they are in the file (CN is NA where it
  # has none). CN_inferred is the copy number implied by the log ratio
  # alone (2 * 2^log.ratio), kept separately: it ignores purity and ploidy,
  # so it isn't a substitute for CN from the callers in use (battenberg,
  # cnvkit, sequenza), but may be useful for a module without CN
  if(!"CN" %in% colnames(all_segs)) all_segs$CN = NA_real_
  all_segs$CN_inferred = 2 * 2^all_segs$log.ratio
  if(max_CN > 0){
    all_segs = dplyr::mutate(all_segs, CN = ifelse(CN > max_CN, max_CN, CN),
                             CN_inferred = ifelse(CN_inferred > max_CN, max_CN, CN_inferred))
  }
  # dummy (placeholder) segments are dropped, along with the column, so
  # callers never see them
  if("dummy_segment" %in% colnames(all_segs)){
    if(verbose){
      message(sprintf("dropping %d dummy segment row(s)",
                      sum(all_segs$dummy_segment %in% 1)))
    }
    all_segs = dplyr::filter(all_segs, !dummy_segment %in% 1) %>%
      dplyr::select(-dummy_segment)
  }
  no_cn = is.na(all_segs$CN)
  if(any(no_cn)){
    by_module = if("module" %in% colnames(all_segs)){
      counts = table(all_segs$module[no_cn], useNA = "ifany")
      paste0(" (module: ", paste(names(counts), counts, sep = " ", collapse = ", "), ")")
    }else ""
    message(sprintf("%d segment(s) in %d sample(s) have no CN%s; CN_inferred has 2 * 2^log.ratio for them",
                    sum(no_cn), dplyr::n_distinct(all_segs$ID[no_cn]), by_module))
  }
  # segments of a sample (profile) should never overlap
  all_segs <- GAMBLR.utils::handle_segment_overlaps(all_segs, action = overlap_action)
  # return S3 class with CN segments and genome_build
  all_segs <- create_seg_data(all_segs, projection)
  return(all_segs)
}

# Column types for a combined seg file, by name, so they're right whatever
# the column order (module, the pipeline each segment came from, is text)
.combined_seg_col_types = function(path){
  header = colnames(suppressMessages(readr::read_tsv(path, n_max = 0, progress = FALSE)))
  types = c(ID = "c", chrom = "c", start = "i", end = "i", module = "c", LOH_flag = "i",
            CN = "d", log.ratio = "d", seg.mean = "d", dummy_segment = "i")
  types = types[names(types) %in% header]
  do.call(readr::cols, c(as.list(types), list(.default = readr::col_guess())))
}
