#' @title Get each sample's segment coverage and baseline copy number
#'
#' @description For the samples in `these_samples_metadata` (or every sample
#' in a seg file), return how much of the autosomal genome their copy number
#' segments cover and their baseline copy number, from a small cache built
#' once per seg file. Use it to find samples whose segments are too sparse
#' for a reliable baseline (`ok` is FALSE), e.g. to remove them before
#' [get_cn_segments] with `baseline_method` or [get_data_by_region].
#'
#' @details The cache is a tab-separated file named after the seg file
#' (`<seg file>.baselines.tsv`): next to the seg file when that directory is
#' writable (shared by everyone using it), otherwise in the per-user cache
#' directory (`tools::R_user_dir("GAMBLR", "cache")`). It is (re)built when
#' it is missing or when the seg file was modified after it, from a read of
#' only the columns needed (ID, chrom, start, end, log.ratio or CN, and
#' dummy_segment). Placeholder (dummy) segments are left out, as in
#' [get_cn_segments]; copy number is the file's CN column (segments
#' without CN are left out, never backfilled from log.ratio), counted up to 20
#' and over autosomes only, as for every GAMBLR baseline (so these match
#' [GAMBLR.utils::segmented_data_to_ploidy_estimate] and
#' [GAMBLR.utils::add_segment_baseline] for the same segments).
#' Coverage and baselines are as in [GAMBLR.utils::check_segment_coverage]
#' and [GAMBLR.utils::add_segment_baseline] (methods "mode",
#' "mode_rounded" and "mean"; see [GAMBLR.utils::baseline_methods]).
#'
#' @param these_samples_metadata Metadata (sample_id, seq_type): each sample
#' is looked up in the combined seg file for its seq_type (genome or
#' capture), or in the per-pipeline seg files with `pipeline`. Give this or
#' `seg_file`.
#' @param projection Genome build ("grch37" or "hg38").
#' @param seg_file Alternatively, the path to a seg file: all its samples.
#' @param pipeline,priority Optional: look samples up in the per-pipeline seg
#' files instead of the combined one, as in [get_cn_segments]. Each sample is
#' taken from the first pipeline (in `priority`, then `pipeline` order) whose
#' file has it, and the result has a `pipeline` column naming it.
#' @param min_genome_coverage,min_autosomes Thresholds for `ok`, as used by
#' [get_cn_segments]. Defaults 0.5 and 20.
#' @param refresh Set to TRUE to rebuild the cache(s).
#' @param verbose Report where caches are read from or written to.
#'
#' @return A data frame with one row per sample: sample_id, seq_type (with
#' metadata), has_segments, coverage, n_autosomes, n_beyond_chrom_end
#' (segments past the chromosome ends of `projection`, e.g. coordinates
#' from another build), n_no_cn (segments without CN), fixed_baseline
#' (TRUE for controlfreec samples, whose baselines are 2 and which don't
#' need to cover the genome; see [GAMBLR.utils::add_segment_baseline]),
#' baseline_mode,
#' baseline_mode_rounded, baseline_mean, coverage_ok, build_ok, cn_ok (the
#' segments have CN, so there is a baseline) and ok (all three; NA for
#' samples without segments). With `pipeline`, also pipeline (NA for samples
#' in none of the files).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' b = get_cn_baselines(meta)
#' too_sparse = b$sample_id[b$ok %in% FALSE]
#' meta_ok = dplyr::filter(meta, !sample_id %in% too_sparse)
#' }
get_cn_baselines = function(these_samples_metadata = NULL, projection = "grch37",
                            seg_file = NULL, pipeline = NULL, priority = NULL,
                            min_genome_coverage = 0.5, min_autosomes = 20,
                            refresh = FALSE, verbose = FALSE){
  flag = function(b){
    # NA for samples without segments (nothing to judge)
    # samples with a fixed baseline (controlfreec) don't need to cover the genome
    b$coverage_ok = ifelse(b$has_segments, b$fixed_baseline %in% TRUE |
                             (b$coverage >= min_genome_coverage & b$n_autosomes >= min_autosomes), NA)
    b$build_ok = ifelse(b$has_segments, b$n_beyond_chrom_end == 0, NA)
    # a baseline needs CN (no backfilling from log.ratio)
    b$cn_ok = ifelse(b$has_segments, !is.na(b$baseline_mode), NA)
    b$ok = b$coverage_ok & b$build_ok & b$cn_ok
    b
  }
  if(!is.null(seg_file)){
    b = .cn_baselines_for_file(seg_file, projection, refresh, verbose)
    b = data.frame(sample_id = b$ID, has_segments = TRUE, b[, -1, drop = FALSE])
    return(flag(b))
  }
  if(is.null(these_samples_metadata)) stop("give these_samples_metadata or seg_file")
  meta = as.data.frame(these_samples_metadata)
  meta = meta[meta$seq_type %in% c("genome", "capture"), c("sample_id", "seq_type"), drop = FALSE]
  order = if(is.null(pipeline)) NULL else .pipeline_order(pipeline, priority)
  out = do.call(rbind, lapply(unique(meta$seq_type), function(st){
    m = meta[meta$seq_type == st, , drop = FALSE]
    if(is.null(order)){
      b = .cn_baselines_for_file(.cnv_combined_path(st, projection), projection,
                                 refresh, verbose)
      from = NULL
    }else{
      # each sample from the first pipeline whose file has it
      b = NULL
      from = rep(NA_character_, nrow(m))
      for(p in order){
        path = .cnv_pipeline_path(p, st, projection)
        if(!file.exists(path)) next
        bp = .cn_baselines_for_file(path, projection, refresh, verbose)
        take = is.na(from) & m$sample_id %in% bp$ID
        from[take] = p
        b = rbind(b, bp[bp$ID %in% m$sample_id[take], , drop = FALSE])
      }
      if(is.null(b)) b = data.frame(ID = character(0))
    }
    at = match(m$sample_id, b$ID)
    out_st = data.frame(sample_id = m$sample_id, seq_type = st, has_segments = !is.na(at),
               coverage = b$coverage[at], n_autosomes = b$n_autosomes[at],
               n_beyond_chrom_end = b$n_beyond_chrom_end[at], n_no_cn = b$n_no_cn[at],
               fixed_baseline = b$fixed_baseline[at],
               baseline_mode = b$baseline_mode[at], baseline_mode_rounded = b$baseline_mode_rounded[at],
               baseline_mean = b$baseline_mean[at],
               stringsAsFactors = FALSE)
    if(!is.null(from)) out_st$pipeline = from
    out_st
  }))
  rownames(out) = NULL
  flag(out)
}

# Path of a per-pipeline seg file (written by the GAMBL per-pipeline merge
# script): the per_pipeline directory next to the combined merge in
# results_merged$cnv
.cnv_pipeline_path = function(pipeline, seq_type, projection){
  merges = dirname(check_config_and_value("results_merged$cnv"))
  paste0(check_config_and_value("project_base"), merges, "/per_pipeline/",
         pipeline, "--", seq_type, "--", projection, ".seg")
}

# Pipelines in the order they are used: those in priority first (in that
# order), then the rest of pipeline
.pipeline_order = function(pipeline, priority = NULL){
  if(length(pipeline) == 0) stop("pipeline must name at least one pipeline")
  bad = setdiff(priority, pipeline)
  if(length(bad) > 0){
    stop("priority names pipelines not in pipeline: ", paste(bad, collapse = ", "))
  }
  unique(c(priority, pipeline))
}

# Path of the combined seg file for a seq_type and build (icgc_dart merge,
# or the gambl one without read access), as in get_cn_segments
.cnv_combined_path = function(seq_type, projection){
  template = check_config_and_value("results_flatfiles$cnv_combined$icgc_dart")
  path = paste0(check_config_and_value("project_base"), glue::glue(template))
  if(file.access(path, 4) == -1){
    template = check_config_and_value("results_flatfiles$cnv_combined$gambl")
    path = paste0(check_config_and_value("project_base"), glue::glue(template))
  }
  path
}

.cn_baselines_for_file = function(seg_file, projection = "grch37",
                                  refresh = FALSE, verbose = FALSE){
  if(!file.exists(seg_file)) stop("seg file not found: ", seg_file)
  shared = paste0(seg_file, ".baselines.tsv")
  user_dir = tools::R_user_dir("GAMBLR", "cache")
  user = file.path(user_dir, paste0(gsub("[^A-Za-z0-9._-]", "_", seg_file), ".baselines.tsv"))
  fresh = function(f) file.exists(f) && file.mtime(f) >= file.mtime(seg_file)
  if(!refresh){
    for(f in c(shared, user)){
      if(fresh(f)){
        cached = as.data.frame(readr::read_tsv(f, col_types = readr::cols(ID = "c", fixed_baseline = "l",
                                                                          .default = "d"),
                                               progress = FALSE))
        # caches from before the build check or the three baseline methods
        # were added are rebuilt
        if(all(c("n_beyond_chrom_end", "baseline_mode_rounded", "n_no_cn", "fixed_baseline") %in%
               colnames(cached))){
          if(verbose) message("baselines from cache: ", f)
          return(cached)
        }
      }
    }
  }
  if(verbose) message("computing baselines for ", basename(seg_file), " (cached for later use)")
  table = .compute_cn_baselines(seg_file, projection)
  target = if(file.access(dirname(seg_file), 2) == 0) shared else {
    dir.create(user_dir, recursive = TRUE, showWarnings = FALSE); user
  }
  ok = tryCatch({ readr::write_tsv(table, target); TRUE }, error = function(e) FALSE)
  if(verbose) message(if(ok) paste("baselines cached in", target) else "baselines could not be cached")
  table
}

# Coverage and baselines for every sample of a seg file, reading only the
# columns needed
.compute_cn_baselines = function(seg_file, projection){
  header = colnames(suppressMessages(readr::read_tsv(seg_file, n_max = 0, progress = FALSE)))
  cols = intersect(c("ID", "chrom", "start", "end", "module", "log.ratio", "seg.mean", "CN",
                    "dummy_segment"), header)
  s = suppressMessages(readr::read_tsv(seg_file, col_select = dplyr::all_of(cols),
                                       na = c("NA", "NaN"), lazy = TRUE, progress = FALSE))
  s = as.data.frame(s)
  if(!"log.ratio" %in% colnames(s) && "seg.mean" %in% colnames(s)) s$log.ratio = s$seg.mean
  if("dummy_segment" %in% colnames(s)){
    s = s[is.na(s$dummy_segment) | s$dummy_segment != 1, , drop = FALSE]
  }
  # the file's CN only (as get_cn_segments returns it); segments without
  # CN don't count, and samples without any have no baseline
  cn = if("CN" %in% colnames(s)) as.numeric(s$CN) else rep(NA_real_, nrow(s))
  n_no_cn = table(factor(s$ID[is.na(cn)], levels = unique(s$ID)))
  # copy number is capped (at 20) by the shared estimator
  auto = gsub("^chr", "", as.character(s$chrom)) %in% as.character(1:22) & !is.na(cn)
  cov = GAMBLR.utils::check_segment_coverage(s[, c("ID", "chrom", "start", "end")],
                                             genome_build = projection)
  len = as.numeric(s$end[auto]) - as.numeric(s$start[auto])
  id = as.character(s$ID[auto])
  mode = GAMBLR.utils:::.baseline_from(id, cn[auto], len, "mode")
  mode_rounded = GAMBLR.utils:::.baseline_from(id, cn[auto], len, "mode_rounded")
  mean = GAMBLR.utils:::.baseline_from(id, cn[auto], len, "mean")
  # samples from modules with a fixed baseline (controlfreec: 2)
  fixed = GAMBLR.utils:::.fixed_baselines(s)
  mode[names(fixed)] = fixed; mode_rounded[names(fixed)] = fixed; mean[names(fixed)] = fixed
  # segments past the chromosome ends of the build (e.g. coordinates from
  # another build), as in GAMBLR.utils::check_segments
  chrom_len = GAMBLR.utils:::.chrom_lengths(projection)[gsub("^chr", "", as.character(s$chrom))]
  beyond = !is.na(chrom_len) & s$end > chrom_len + 1e5
  n_beyond = table(factor(s$ID[beyond], levels = cov$ID))
  data.frame(ID = cov$ID, coverage = round(cov$coverage, 4), n_autosomes = cov$n_autosomes,
             n_beyond_chrom_end = as.integer(n_beyond),
             n_no_cn = as.integer(n_no_cn[cov$ID]),
             fixed_baseline = cov$ID %in% names(fixed),
             baseline_mode = round(as.numeric(mode[cov$ID]), 4),
             baseline_mode_rounded = as.numeric(mode_rounded[cov$ID]),
             baseline_mean = round(as.numeric(mean[cov$ID]), 4),
             stringsAsFactors = FALSE)
}
