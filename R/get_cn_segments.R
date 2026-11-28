


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
#' @param pipeline Optional: get segments from specific pipelines instead of
#' the combined merge (one chosen source per sample). Names as in the
#' per-pipeline merge files (the `per_pipeline` directory next to the
#' combined merge in `results_merged$cnv`, named
#' `<pipeline>--<seq_type>--<projection>.seg`), e.g. "battenberg",
#' "controlfreec", "sequenza", "cnvkit", "purecn" or "gatk_cnv" (GATK
#' ModelSegments; capture only, with CN = round(2 * 2^log.ratio), i.e. not
#' purity-corrected). Published segmentations bundled in GAMBLR.data can be
#' used the same way, for comparison with GAMBL's: "dlbclass" (DLBClass,
#' Chapuy et al.; CN is ABSOLUTE's purity/ploidy-corrected copy number) and
#' "schmitz_array" (Schmitz et al. arrays), for capture samples of
#' dlbcl_chapuy and dlbcl_schmitz (see GAMBLR.data::published_cn_sources).
#' "all" uses every pipeline with a per-pipeline file, and the published
#' sources
#' for the seq_types in the metadata (see [cn_pipeline_availability]), and
#' then `one_per_sample` defaults to FALSE. Default NULL uses the combined
#' merge. A pipeline without a file for a seq_type is skipped, with a
#' message (silently with "all").
#' @param priority With several pipelines, the ones to prefer, in order
#' (with `source_rule = "preference"`, the order of preference);
#' pipelines not named follow in `pipeline` order. E.g. `pipeline =
#' c("purecn", "battenberg"), priority = "battenberg"` takes each sample from
#' battenberg when it has it, and from purecn otherwise. Default NULL: the
#' order of `pipeline`.
#' @param one_per_sample With several pipelines, TRUE (default, except with
#' `pipeline = "all"`) returns each sample's segments from one pipeline only
#' (the first, by `priority`, whose file has the sample); FALSE returns every
#' pipeline's segments for every sample, e.g. to compare them (the overlap
#' check is then per pipeline).
#' @param source_rule How each sample's pipeline is chosen:
#' \itemize{
#'   \item "cnv_master" (default): the combined merge, whose pipeline per
#'     sample is set by a fixed priority when the merges are built.
#'   \item "qc_score": each sample's lowest-scoring GAMBL pipeline by
#'     [GAMBLR.utils::compare_segmentations]. Every sample in the
#'     per-pipeline files of the seq_type and `projection` is scored, not
#'     only those in `these_samples_metadata` (which only subsets the result
#'     afterwards), so the choice depends on the files alone and is the same
#'     for every caller. A segmentation that can't give a genome-wide
#'     baseline (see [get_cn_baselines]) is never chosen, so samples whose
#'     segmentations all fail that are dropped. controlfreec, run on genomes
#'     without a matched normal, reports only altered segments: its scores
#'     carry a large coverage penalty and aren't comparable with other
#'     pipelines', but it is usually a sample's only candidate. The scores are cached in a table next to the files
#'     (`cn_qc_scores--<seq_type>--<projection>.tsv`), or in the per-user cache
#'     directory, rebuilt when a file is newer, the pipelines available
#'     change, or GAMBLR.utils's version changes. A choice can change when a
#'     merge is regenerated or a pipeline is added: to keep it, save
#'     `attr(result, "selection")` and pass it back as `source_rule`.
#'   \item "preference": each sample's first pipeline, in order of
#'     preference (`priority`; by default c("purecn", "gatk_cnv", "cnvkit")
#'     for capture and c("battenberg", "sequenza", "controlfreec") for
#'     genomes, purity-corrected pipelines first; pipelines not named follow,
#'     alphabetically), that passes the QC gates of the scores used by
#'     "qc_score": it can give a genome-wide baseline, has no contradictory
#'     deep loss, and isn't over-segmented (more segments per 100 Mb than
#'     `max_oversegmentation` times the median of the same pipeline, so a
#'     pipeline that segments finely by design isn't held to another's
#'     density). A sample whose pipelines all fail a gate gets the eligible
#'     one closest to passing (no contradiction first, then the least
#'     over-segmented for its pipeline). The selection records each sample's skipped pipelines
#'     and why (`skipped`) and whether its pipeline passed every gate
#'     (`all_gates_passed`). On LymphGen's training data (Schmitz exomes),
#'     this beats "qc_score", which prefers gatk_cnv for almost every sample
#'     because PureCN segments more finely.
#'   \item A saved selection: a table (e.g. `attr(result, "selection")`, or
#'     [GAMBLR.utils::compare_segmentations]'s selected rows) or the path of
#'     a TSV of one, with sample_id (or ID), pipeline and optionally seq_type.
#' }
#' With anything but "cnv_master", samples without a chosen pipeline are
#' dropped, with a message; `pipeline` can't be given too.
#' @param max_oversegmentation With `source_rule = "preference"`, the
#' most segments per 100 Mb a pipeline may have, as a multiple of that
#' pipeline's median (per seq_type and native genome version), before the
#' next preferred pipeline is used. Default 2; NA turns the gate off.
#' @param calibrate Convert the segments of samples whose pipeline has
#' relative copy number (e.g. gatk_cnv and CNVkit, whose CN is only a
#' rescaled log ratio; see [GAMBLR.utils::cn_pipeline_types]) to absolute
#' copy number with [GAMBLR.utils::calibrate_relative_cn], using the same
#' sample's segmentation from the first pipeline in `calibrate_reference`
#' that passes QC (eligible in the QC scores of `source_rule = "qc_score"`,
#' i.e. it gives a baseline, with no contradiction penalty, i.e. not a failed
#' absolute fit). The fit uses genome-wide segments, also with `region`.
#' Calibrated segments get `cn_type` "absolute" (the original CN is in
#' `CN_relative`), so copy number features are then called on CN (e.g. by
#' GAMBLR.predict::lymphgen_input); samples without a passing reference, or
#' whose fit is flagged, stay relative. Each sample's reference and fit are
#' in `attr(result, "calibration")`. Default FALSE.
#' @param calibrate_reference Pipelines with absolute copy number to
#' calibrate against, in order of preference. Default c("purecn",
#' "battenberg").
#' @param flavour Deprecated: use `pipeline`. "combined" (default) is the
#' combined merge; any other value is used as `pipeline`.
#' @param max_CN Cap CN values at this value. Default is 20. This is to avoid extreme CN values dominating the colour scale in visualisations.
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
#' with [annotate_cn_availability] and `cn_available`. With `pipeline`, each
#' sample's baseline comes from the pipeline its segments come from.
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
#' have usable copy number data (and, with `pipeline`, `cn_pipeline`: the
#' pipeline each sample's data come from, or with `one_per_sample = FALSE`
#' a comma-separated list of the pipelines with segments for it). See
#' [annotate_data_availability].
#' With the combined merge (`source_rule = "cnv_master"`), `cn_pipeline` is
#' the module each sample's segments come from (NA without segments).
#' With a `source_rule` other than "cnv_master", `cn_pipeline` is each
#' sample's chosen pipeline (and `cn_qc_score` its score, when known).
#' Default FALSE.
#' @return A data frame with CN segments for the specified region. Dummy
#' (placeholder) segments are dropped, and the `dummy_segment` column with
#' them. Includes `module` (the pipeline each segment came from, e.g.
#' battenberg, cnvkit, sequenza) when the seg file has it, and a `pipeline`
#' column naming each segment's source: the per-pipeline merge it was read
#' from, or for the combined merge a copy of `module`. With a `source_rule`
#' other than "cnv_master", each sample's chosen pipeline (and its scores)
#' are in `attr(result, "selection")`. `CN` and
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
#' \dontrun{
#' # segments from one pipeline only
#' purecn_segs = get_cn_segments(capture_metadata, pipeline = "purecn")
#'
#' # battenberg where available, purecn otherwise
#' segs = get_cn_segments(mixed_seq_type_meta, pipeline = c("purecn", "battenberg"),
#'                        priority = "battenberg")
#'
#' # every pipeline's segments, side by side
#' all_segs = get_cn_segments(capture_metadata, pipeline = c("cnvkit", "purecn"),
#'                            one_per_sample = FALSE)
#'
#' # every available pipeline, side by side
#' all_segs = get_cn_segments(capture_metadata, pipeline = "all")
#'
#' # each sample's pipeline chosen by QC score (scored over all samples in
#' # the files); save the choice to reuse it later
#' qc_segs = get_cn_segments(capture_metadata, source_rule = "qc_score")
#' readr::write_tsv(attr(qc_segs, "selection"), "cn_selection.tsv")
#' same_segs = get_cn_segments(capture_metadata, source_rule = "cn_selection.tsv")
#' }
#'

get_cn_segments <- function(these_samples_metadata,
                            projection = "grch37",
                            pipeline = NULL,
                            priority = NULL,
                            one_per_sample = TRUE,
                            source_rule = "cnv_master",
                            max_oversegmentation = 2,
                            calibrate = FALSE,
                            calibrate_reference = c("purecn", "battenberg"),
                            flavour = "combined",
                            this_seq_type,
                            max_CN = 20,
                            verbose = FALSE,
                            overlap_action = "warn",
                            region = NULL,
                            baseline_method = "mode",
                            with_metadata = FALSE) {
  # flavour is deprecated: anything but "combined" is a pipeline
  if(flavour != "combined"){
    if(!is.null(pipeline)) stop("give pipeline or flavour, not both (flavour is deprecated)")
    message("flavour is deprecated; use pipeline = \"", flavour, "\"")
    pipeline <- flavour
  }
  # a source_rule other than cnv_master: each sample's chosen pipeline, then
  # the chosen pipelines read side by side and each sample's kept (below)
  chosen <- NULL
  if(!identical(source_rule, "cnv_master")){
    if(!is.null(pipeline)) stop("give pipeline or source_rule, not both")
    st <- if(missing(these_samples_metadata) || is.null(these_samples_metadata)) c("genome", "capture") else
      intersect(unique(these_samples_metadata$seq_type), c("genome", "capture"))
    chosen <- .cn_source_selection(source_rule, st, projection, verbose, priority = priority,
                                   max_oversegmentation = max_oversegmentation)
    if(!missing(these_samples_metadata) && !is.null(these_samples_metadata)){
      mine <- .selection_for(chosen, these_samples_metadata$sample_id, these_samples_metadata$seq_type)
      chosen <- chosen[paste(chosen$sample_id, chosen$seq_type) %in% paste(mine$sample_id, mine$seq_type), , drop = FALSE]
    }
    if(nrow(chosen) == 0) stop("source_rule chose no pipeline for any of these samples")
    pipeline <- unique(chosen$pipeline)
    one_per_sample <- FALSE
  }
  quiet_missing <- identical(pipeline, "all") || !is.null(chosen)
  if(identical(pipeline, "all")){
    # every pipeline with a per-pipeline file for these seq_types; compared
    # side by side unless one_per_sample is given
    st <- if(missing(these_samples_metadata) || is.null(these_samples_metadata)) c("genome", "capture") else
      intersect(unique(these_samples_metadata$seq_type), c("genome", "capture"))
    r <- .resolve_cn_pipelines(pipeline, priority, st, projection)
    pipeline <- r$pipeline; priority <- r$priority
    if(missing(one_per_sample)) one_per_sample <- FALSE
  }
  if(!is.null(pipeline)) pipeline_order <- .pipeline_order(pipeline, priority)
  if(with_metadata){
    if(missing(these_samples_metadata) || is.null(these_samples_metadata)){
      stop("with_metadata = TRUE needs these_samples_metadata")
    }
    call <- match.call()
    if(is.null(chosen)){
      call$pipeline <- pipeline; call$priority <- priority; call$one_per_sample <- one_per_sample
    }
    return(.with_metadata(call, parent.frame(), "seg_data", function(){
      # with a source_rule, cn_pipeline is each sample's chosen pipeline, and
      # cn_available whether that pipeline's segments give a baseline
      if(!is.null(chosen)){
        a <- these_samples_metadata
        own <- .selection_for(chosen, a$sample_id, a$seq_type)
        a$cn_pipeline <- own$pipeline
        if("score" %in% colnames(own)) a$cn_qc_score <- own$score
        ok <- rep(FALSE, nrow(a))
        for(p in unique(stats::na.omit(a$cn_pipeline))){
          rows <- which(a$cn_pipeline %in% p)
          b <- get_cn_baselines(a[rows, , drop = FALSE], projection = projection, pipeline = p)
          ok[rows] <- b$ok[match(paste(a$sample_id[rows], a$seq_type[rows]),
                                 paste(b$sample_id, b$seq_type))] %in% TRUE
        }
        return(.set_available(a, "cn_available", ok, "copy number data", FALSE, TRUE))
      }
      a <- annotate_cn_availability(these_samples_metadata, projection = projection,
                                    pipeline = pipeline, priority = priority)
      # side by side, cn_pipeline lists every pipeline with segments for the sample
      if(!is.null(pipeline) && !one_per_sample){
        av <- cn_pipeline_availability(these_samples_metadata, projection = projection,
                                       pipeline = pipeline)
        has <- as.matrix(av[, pipeline, drop = FALSE]) > 0
        listed <- apply(has, 1, function(x) paste(pipeline[x %in% TRUE], collapse = ","))
        a$cn_pipeline <- listed[match(paste(a$sample_id, a$seq_type), paste(av$sample_id, av$seq_type))]
        a$cn_pipeline[a$cn_pipeline %in% ""] <- NA
      }
      # the combined merge: cn_pipeline is the module cnv_master took each
      # sample's segments from
      if(is.null(pipeline)) a$cn_pipeline <- .combined_modules(a, projection)
      a
    }))
  }
  if (!missing(this_seq_type)) {
    stop("this_seq_type is deprecated. Subset your metadata instead.")
  }
  # optional baselines: each sample's genome-wide baseline copy number, from
  # the cached per-file table (see get_cn_baselines)
  add_baselines = function(x, path, table = NULL) x
  # samples without a baseline, from every file read, reported once at the end
  no_baseline <- new.env()
  no_baseline$ids <- character(0)
  # baselines come from each seg file's cache (combined or per-pipeline)
  if (!is.null(baseline_method)) {
    baseline_method <- match.arg(baseline_method, GAMBLR.utils::baseline_methods())
    # baselines from a seg file's cache, or a given table (published sources);
    # path labels the messages
    add_baselines = function(x, path, table = NULL) {
      b <- if(is.null(table)) .cn_baselines_for_file(path, projection = projection, verbose = verbose) else table
      # samples that can't get a reliable baseline (as in get_cn_baselines:
      # too sparse, coordinates past the chromosome ends, or no CN) get NA
      # rather than stopping the call; functions that need a baseline stop
      # on them, naming them
      b <- b[b$ID %in% unique(x$ID), , drop = FALSE]
      bad <- (!b$fixed_baseline %in% TRUE & (b$coverage < 0.5 | b$n_autosomes < 20)) |
        b$n_beyond_chrom_end > 0 | is.na(b$baseline_mode)
      b[bad, paste0("baseline_", GAMBLR.utils::baseline_methods())] <- NA
      missing_ids <- setdiff(unique(x$ID), b$ID)
      if (verbose && any(bad)) {
        message(basename(path), ": ", sum(bad), " sample(s) have no baseline_cn (segments too sparse, past the ",
                "chromosome ends, or without CN): ", paste(b$ID[bad], collapse = ", "))
      }
      if (verbose && length(missing_ids) > 0) {
        message(basename(path), ": ", length(missing_ids), " sample(s) are missing from the baseline cache: ",
                paste(missing_ids, collapse = ", "))
      }
      no_baseline$ids <- c(no_baseline$ids, b$ID[bad], missing_ids)
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
  if(!is.null(chosen)){
    n_none <- 0
    if(!missing(these_samples_metadata)){
      these_samples_metadata <- dplyr::filter(these_samples_metadata, seq_type %in% c("genome", "capture"))
      own <- .selection_for(chosen, these_samples_metadata$sample_id, these_samples_metadata$seq_type)
      none <- is.na(own$pipeline)
      n_none <- length(unique(these_samples_metadata$sample_id[none]))
      if(verbose && n_none > 0){
        message("no pipeline chosen for: ", paste(unique(these_samples_metadata$sample_id[none]), collapse = ", "))
      }
      these_samples_metadata <- these_samples_metadata[!none, , drop = FALSE]
    }
    counts <- sort(table(chosen$pipeline), decreasing = TRUE)
    message("pipelines ", if(identical(source_rule, "qc_score")) "by QC score" else
              if(identical(source_rule, "preference")) "by preference and QC gates" else "from the saved selection",
            " for ", nrow(chosen), " sample(s): ", paste(names(counts), counts, collapse = ", "),
            if(n_none > 0) paste0("; ", n_none, " sample(s) with none are dropped") else "")
    if(identical(source_rule, "preference") && any(!chosen$all_gates_passed)){
      message(sum(!chosen$all_gates_passed), " sample(s) had no pipeline passing every QC gate; the eligible ",
              "pipeline closest to passing is used (see attr(<result>, \"selection\"))")
    }
    if(is.character(source_rule) && length(source_rule) == 1 && source_rule %in% c("qc_score", "preference") &&
       !isTRUE(getOption("GAMBLR.results.shown_selection_tip"))){
      message("(the choice can change when the merges change; to keep it, save attr(<result>, \"selection\") ",
              "and pass it as source_rule)")
      options(GAMBLR.results.shown_selection_tip = TRUE)
    }
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
  if (!is.null(chosen)) {
    # each pipeline read only for the samples it was chosen for
    parts <- lapply(pipeline, function(p){
      mine <- chosen[chosen$pipeline == p, , drop = FALSE]
      ids <- list(genome = intersect(genome_ids, mine$sample_id[mine$seq_type %in% c("genome", NA)]),
                  capture = intersect(capture_ids, mine$sample_id[mine$seq_type %in% c("capture", NA)]))
      tryCatch(.get_pipeline_segments(p, FALSE, projection, ids, in_region, add_baselines, verbose,
                                      quiet_missing = TRUE, summarise = verbose),
               error = function(e) if(grepl("^no segments found", conditionMessage(e))) NULL else stop(e))
    })
    all_segs <- dplyr::bind_rows(parts)
    if(nrow(all_segs) == 0) stop("no segments found for these samples in their chosen pipelines")
    absent <- setdiff(these_samples_metadata$sample_id, all_segs$ID)
    if(length(absent)){
      message(length(absent), " sample(s) have no segments from their chosen pipeline",
              if(verbose) paste0(": ", paste(absent, collapse = ", ")) else
                if(!is.null(region)) " in the region" else "")
    }
  } else if (!is.null(pipeline)) {
    all_segs <- .get_pipeline_segments(pipeline_order, one_per_sample, projection,
                                       list(genome = genome_ids, capture = capture_ids),
                                       in_region, add_baselines, verbose, quiet_missing)
  } else {
    cnv_flatfile_template <- check_config_and_value(
      "results_flatfiles$cnv_combined$icgc_dart"
    )
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

      # column types from the file's header (see .combined_seg_col_types)
      coltypes <- .combined_seg_col_types(full_cnv_path)
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
        print(nrow(seg))
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
    # the source of each segment, named as in the per-pipeline merges
    if("module" %in% colnames(all_segs) && !"pipeline" %in% colnames(all_segs)) all_segs$pipeline <- all_segs$module
  }
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
  # relative copy number converted to absolute, against a reference
  # segmentation of the same sample
  calibration <- NULL
  if(calibrate){
    cal <- .calibrate_relative_segments(all_segs, calibrate_reference, projection, baseline_method,
                                        max_CN, verbose)
    all_segs <- cal$segs
    calibration <- cal$calibration
  }
  no_cn = is.na(all_segs$CN)
  if(verbose && any(no_cn)){
    by_module = if("module" %in% colnames(all_segs)){
      counts = table(all_segs$module[no_cn], useNA = "ifany")
      paste0(" (module: ", paste(names(counts), counts, sep = " ", collapse = ", "), ")")
    }else ""
    message(sprintf("%d segment(s) in %d sample(s) have no CN%s; CN_inferred has 2 * 2^log.ratio for them",
                    sum(no_cn), dplyr::n_distinct(all_segs$ID[no_cn]), by_module))
  }
  n_no_baseline <- length(intersect(unique(no_baseline$ids), all_segs$ID))
  if(n_no_baseline > 0){
    message(n_no_baseline, " sample(s) have no baseline_cn, so functions that need one ",
            "(e.g. GAMBLR.utils::center_segments) stop on them; with_metadata = TRUE marks them ",
            "(cn_available)", if(!verbose) "; verbose = TRUE lists them" else "")
  }
  # segments of a sample (profile) should never overlap
  # (with several pipelines kept side by side, per pipeline)
  overlap_by <- if (!is.null(pipeline) && !one_per_sample) {
    c("ID", intersect("seg_seq_type", colnames(all_segs)), "pipeline")
  } else NULL
  all_segs <- GAMBLR.utils::handle_segment_overlaps(all_segs, action = overlap_action,
                                                    by = overlap_by, list_samples = verbose)
  # return S3 class with CN segments and genome_build
  all_segs <- create_seg_data(all_segs, projection)
  if(!is.null(chosen)){
    attr(all_segs, "selection") <- chosen[chosen$sample_id %in% unique(all_segs$ID), , drop = FALSE]
  }
  if(!is.null(calibration)) attr(all_segs, "calibration") <- calibration
  return(all_segs)
}

# The module (pipeline) the combined merge takes each sample's segments from,
# read from the ID and module columns of its file for each seq_type (several
# are comma-separated; NA without segments)
.combined_modules = function(meta, projection){
  out = rep(NA_character_, nrow(meta))
  for(st in intersect(unique(meta$seq_type), c("genome", "capture"))){
    path = .cnv_combined_path(st, projection)
    if(!file.exists(path)) next
    header = colnames(suppressMessages(readr::read_tsv(path, n_max = 0, progress = FALSE)))
    if(!"module" %in% header) next
    m = unique(as.data.frame(suppressMessages(readr::read_tsv(path, col_select = c("ID", "module"),
                                                              col_types = "cc", progress = FALSE))))
    m = m[!is.na(m$module), , drop = FALSE]
    per_id = tapply(m$module, m$ID, function(x) paste(sort(unique(x)), collapse = ","))
    rows = which(meta$seq_type == st)
    out[rows] = as.character(per_id[meta$sample_id[rows]])
  }
  out
}

# Segments from the per-pipeline merges, for get_cn_segments(pipeline = ...).
# For each seq_type, the pipelines' files are read in `order`; with
# one_per_sample, each sample comes from the first file that has it (judged
# on the whole file, before any region filter, so a sample isn't switched to
# another pipeline just for lacking segments in the region)
.get_pipeline_segments = function(order, one_per_sample, projection, ids_by_seq_type,
                                  in_region, add_baselines, verbose, quiet_missing = FALSE,
                                  summarise = TRUE){
  out = list()
  for(seq_type in names(ids_by_seq_type)){
    ids = ids_by_seq_type[[seq_type]]
    if(length(ids) == 0) next
    taken = character(0)
    counts = c()
    for(p in order){
      published = p %in% .published_cn_sources()
      path = if(published) p else .cnv_pipeline_path(p, seq_type, projection)
      if(!.pipeline_has_data(p, seq_type, projection)){
        if(!quiet_missing) message("no ", p, " segments for ", seq_type, " samples",
                                   if(published) " (published sources cover capture samples)" else
                                     paste0(" (", path, " not found)"))
        next
      }
      if(verbose) print(path)
      seg = if(published) .published_segments(p, projection, sample_ids = ids) else
        suppressMessages(readr::read_tsv(path, col_types = .combined_seg_col_types(path),
                                         na = c("NA", "NaN"), progress = FALSE, lazy = TRUE))
      seg = seg[seg$ID %in% ids, , drop = FALSE]
      if(one_per_sample){
        seg = seg[!seg$ID %in% taken, , drop = FALSE]
        taken = c(taken, unique(seg$ID))
      }
      counts[p] = dplyr::n_distinct(seg$ID)
      if(nrow(seg) == 0) next
      seg = add_baselines(in_region(seg), path,
                          if(published) .pipeline_baselines(p, seq_type, projection) else NULL)
      seg$pipeline = rep(p, nrow(seg))
      seg$seg_seq_type = rep(seq_type, nrow(seg))
      out[[paste(seq_type, p)]] = seg
    }
    if(summarise && length(counts) > 0){
      found = if(one_per_sample) length(taken) else
        dplyr::n_distinct(unlist(lapply(out[grepl(paste0("^", seq_type, " "), names(out))],
                                        function(x) x$ID)))
      message(seq_type, ": ", paste(counts, "sample(s) from", names(counts), collapse = ", "),
              if(length(ids) > found) paste0("; ", length(ids) - found,
                                               " sample(s) in none of these pipelines") else "")
    }
  }
  if(length(out) == 0){
    stop("no segments found for these samples in pipeline(s): ", paste(order, collapse = ", "))
  }
  dplyr::bind_rows(out)
}

# Column types for a combined seg file, by name, so they're right whatever
# the column order (module, the pipeline each segment came from, is text)
.combined_seg_col_types = function(path){
  header = colnames(suppressMessages(readr::read_tsv(path, n_max = 0, progress = FALSE)))
  # LOH_flag as a number: some modules write it as 0.0/1.0
  types = c(ID = "c", chrom = "c", start = "i", end = "i", module = "c", LOH_flag = "d",
            CN = "d", log.ratio = "d", seg.mean = "d", dummy_segment = "i")
  types = types[names(types) %in% header]
  do.call(readr::cols, c(as.list(types), list(.default = readr::col_guess())))
}
