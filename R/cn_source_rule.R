# source_rule for get_cn_segments: how each sample's pipeline is chosen when
# it isn't the combined merge's (cnv_master) choice

# QC scores (GAMBLR.utils::compare_segmentations) of every sample in every
# GAMBL per-pipeline merge of a seq_type and projection. All samples in the
# files are scored, not just a caller's, so the choice depends on the files
# alone. Cached next to the files (or in the per-user cache directory) and
# rebuilt when a file is newer, the pipelines available change, or
# GAMBLR.utils's version or the scoring rules (.cn_qc_rules) change.
# Published sources are never scored. A segmentation that can't give a
# genome-wide baseline (as in get_cn_baselines) is never selected.
.cn_qc_scores = function(seq_type, projection, verbose = FALSE){
  pipelines = setdiff(.available_cn_pipelines(seq_type, projection), .published_cn_sources())
  if(length(pipelines) == 0) return(NULL)
  files = vapply(pipelines, .cnv_pipeline_path, "", seq_type = seq_type, projection = projection)
  name = paste0("cn_qc_scores--", seq_type, "--", projection, ".tsv")
  shared = file.path(dirname(files[1]), name)
  user_dir = tools::R_user_dir("GAMBLR", "cache")
  user = file.path(user_dir, paste0(gsub("[^A-Za-z0-9._-]", "_", dirname(files[1])), "_", name))
  stamp = c(paste("# pipelines:", paste(pipelines, collapse = ",")),
            paste("# GAMBLR.utils:", as.character(utils::packageVersion("GAMBLR.utils")),
                  "; rules:", .cn_qc_rules))
  newest = max(file.mtime(files))
  col_types = readr::cols(ID = "c", sample_id = "c", seq_type = "c", genome = "c", pipeline = "c",
                          eligible = "l", selected = "l", .default = "d")
  for(f in c(shared, user)){
    if(file.exists(f) && file.mtime(f) >= newest && identical(readLines(f, n = 2), stamp)){
      if(verbose) message("QC scores from cache: ", f)
      return(as.data.frame(readr::read_tsv(f, comment = "#", col_types = col_types, progress = FALSE)))
    }
  }
  message("scoring ", paste(pipelines, collapse = ", "), " for all ", seq_type, " samples (",
          projection, "; first use, then cached)")
  segs = dplyr::bind_rows(lapply(pipelines, function(p){
    path = .cnv_pipeline_path(p, seq_type, projection)
    s = suppressMessages(readr::read_tsv(path, col_types = .combined_seg_col_types(path),
                                         na = c("NA", "NaN"), progress = FALSE))
    if(!"log.ratio" %in% colnames(s) && "seg.mean" %in% colnames(s)) s = dplyr::rename(s, log.ratio = seg.mean)
    if("dummy_segment" %in% colnames(s)) s = s[!s$dummy_segment %in% 1, , drop = FALSE]
    if(!"CN" %in% colnames(s)) s$CN = NA_real_
    s = s[, c("ID", "chrom", "start", "end", "CN", "log.ratio")]
    s$pipeline = rep(p, nrow(s))
    s
  }))
  # strata: the seq_type of the files, and each sample's native genome_build
  # from the metadata (NA, a stratum of its own, for samples it lacks)
  meta = tryCatch(suppressWarnings(suppressMessages(
    get_gambl_metadata(collapse_redundancy = FALSE, include_unavailable = TRUE))),
    error = function(e) NULL)
  ids = unique(segs$ID)
  build = if(is.null(meta)) rep(NA_character_, length(ids)) else {
    m = meta[meta$seq_type == seq_type, , drop = FALSE]
    as.character(m$genome_build[match(ids, m$sample_id)])
  }
  strata = data.frame(sample_id = ids, seq_type = seq_type, genome_build = build, stringsAsFactors = FALSE)
  # eligible: segmentations that can give a baseline, from each file's
  # baseline cache (controlfreec's fixed baseline is exempt from coverage)
  eligible = dplyr::bind_rows(lapply(pipelines, function(p){
    b = .cn_baselines_for_file(.cnv_pipeline_path(p, seq_type, projection), projection, verbose = verbose)
    data.frame(ID = b$ID, pipeline = rep(p, nrow(b)), eligible = .baseline_ok(b), stringsAsFactors = FALSE)
  }))
  scores = GAMBLR.utils::compare_segmentations(segs, strata, genome_build = projection, eligible = eligible)
  scores$seq_type = rep(seq_type, nrow(scores))
  target = if(file.access(dirname(shared), 2) == 0) shared else {
    dir.create(user_dir, recursive = TRUE, showWarnings = FALSE); user
  }
  ok = tryCatch({
    writeLines(stamp, target)
    readr::write_tsv(scores, target, append = TRUE, col_names = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if(verbose) message(if(ok) paste("QC scores cached in", target) else "QC scores could not be cached")
  scores
}

# Scoring rules, part of the cache stamp: change it when the rules change
.cn_qc_rules = "baseline-eligible"

# Whether each sample in a baseline table (from .cn_baselines_for_file) has a
# usable baseline, as in get_cn_baselines: enough of the genome covered
# (unless the module has a fixed baseline), no segments past the chromosome
# ends, and CN
.baseline_ok = function(b, min_genome_coverage = 0.5, min_autosomes = 20){
  (b$fixed_baseline %in% TRUE | (b$coverage >= min_genome_coverage & b$n_autosomes >= min_autosomes)) &
    b$n_beyond_chrom_end == 0 & !is.na(b$baseline_mode)
}

# Each pipeline's place in source_rule = "preference" when none is given:
# purity-corrected pipelines first
.cn_default_preference = list(capture = c("purecn", "gatk_cnv", "cnvkit"),
                              genome = c("battenberg", "sequenza", "controlfreec"))

# source_rule = "preference": each sample's first pipeline, in order of
# preference, that passes the QC gates (from the QC scores): it can give a
# baseline (eligible), has no contradictory deep loss (pen_contradiction),
# and isn't over-segmented, i.e. its segments per 100 Mb are at most
# max_oversegmentation times the median of the same pipeline (per seq_type
# and native genome version), so a pipeline that segments finely by design
# isn't held to another's density. When no pipeline passes every gate, the
# eligible one closest to passing is taken: without a contradiction first,
# then the least over-segmented for its pipeline. Pipelines not in the
# preference follow it, alphabetically.
.cn_preference_selection = function(seq_types, projection, priority = NULL, max_oversegmentation = 2,
                                    verbose = FALSE){
  scores = dplyr::bind_rows(lapply(seq_types, .cn_qc_scores, projection = projection, verbose = verbose))
  if(nrow(scores) == 0) stop("no per-pipeline seg files to score for ", paste(seq_types, collapse = "/"),
                             " in ", projection)
  scores = as.data.frame(scores)
  key = paste(scores$seq_type, scores$genome, scores$pipeline)
  typical = tapply(scores$segments_per_100mb, key, stats::median, na.rm = TRUE)
  scores$oversegmentation = scores$segments_per_100mb / as.numeric(typical[key])
  dense = if(is.null(max_oversegmentation) || is.na(max_oversegmentation)) FALSE else
    !is.na(scores$oversegmentation) & scores$oversegmentation > max_oversegmentation
  scores$gate = ifelse(!scores$eligible %in% TRUE, "no baseline",
                ifelse(!is.na(scores$pen_contradiction) & scores$pen_contradiction > 0, "contradictory deep loss",
                ifelse(dense, "over-segmented", "passed")))
  scores$preference = NA_integer_
  for(st in unique(scores$seq_type)){
    in_st = scores$seq_type == st
    order = if(is.null(priority)) .cn_default_preference[[st]] else priority
    order = c(order, sort(setdiff(unique(scores$pipeline[in_st]), order)))
    scores$preference[in_st] = match(scores$pipeline[in_st], order)
  }
  scores = scores[order(scores$seq_type, scores$sample_id, scores$preference), , drop = FALSE]
  by_sample = split(seq_len(nrow(scores)), paste(scores$seq_type, scores$sample_id))
  picked = vapply(by_sample, function(i){
    ok = i[scores$gate[i] == "passed"]
    if(length(ok) > 0) return(ok[1])
    # none passes: the eligible one closest to passing (no contradiction
    # first, then the least over-segmented for its pipeline)
    eligible = i[scores$eligible[i] %in% TRUE]
    if(length(eligible) == 0) return(NA_integer_)
    contradiction = !is.na(scores$pen_contradiction[eligible]) & scores$pen_contradiction[eligible] > 0
    eligible[order(contradiction, scores$oversegmentation[eligible], scores$preference[eligible])[1]]
  }, 0L)
  # what was passed over for each sample, and why
  skipped = vapply(names(by_sample), function(k){
    i = by_sample[[k]]; j = picked[[k]]
    if(is.na(j)) return(NA_character_)
    before = i[scores$preference[i] < scores$preference[j]]
    if(length(before) == 0) "" else paste0(scores$pipeline[before], " (", scores$gate[before], ")", collapse = "; ")
  }, "")
  sel = scores[stats::na.omit(picked), , drop = FALSE]
  sel$skipped = skipped[!is.na(picked)]
  sel$all_gates_passed = sel$gate == "passed"
  sel$selected = TRUE
  rownames(sel) = NULL
  sel
}

# Each sample's chosen pipeline under a source_rule other than "cnv_master":
# the selected rows of the QC scores ("qc_score"), the first passing
# pipeline in order of preference ("preference"), or a saved selection (a
# table or a path to a TSV of one), as sample_id, seq_type, pipeline (and the
# scores, when known)
.cn_source_selection = function(source_rule, seq_types, projection, verbose = FALSE, priority = NULL,
                                max_oversegmentation = 2){
  if(identical(source_rule, "preference")){
    sel = .cn_preference_selection(seq_types, projection, priority, max_oversegmentation, verbose)
  }else if(identical(source_rule, "qc_score")){
    scores = dplyr::bind_rows(lapply(seq_types, .cn_qc_scores, projection = projection, verbose = verbose))
    if(nrow(scores) == 0) stop("no per-pipeline seg files to score for ", paste(seq_types, collapse = "/"),
                               " in ", projection)
    sel = scores[scores$selected %in% TRUE, , drop = FALSE]
  }else{
    sel = if(is.character(source_rule)){
      if(!file.exists(source_rule)) stop("source_rule must be \"cnv_master\", \"qc_score\", \"preference\", a selection table ",
                                         "or the path of a saved one; not found: ", source_rule)
      readr::read_tsv(source_rule, comment = "#", show_col_types = FALSE, progress = FALSE)
    } else source_rule
    sel = as.data.frame(sel)
    if(!"sample_id" %in% colnames(sel) && "ID" %in% colnames(sel)) sel$sample_id = sel$ID
    missing_cols = setdiff(c("sample_id", "pipeline"), colnames(sel))
    if(length(missing_cols)) stop("source_rule selection lacks: ", paste(missing_cols, collapse = ", "))
    if("selected" %in% colnames(sel)) sel = sel[sel$selected %in% TRUE, , drop = FALSE]
    if(!"seq_type" %in% colnames(sel)) sel$seq_type = NA_character_
  }
  sel$sample_id = as.character(sel$sample_id)
  sel$seq_type = as.character(sel$seq_type)
  sel$pipeline = as.character(sel$pipeline)
  rownames(sel) = NULL
  sel
}

# Rows of a selection for these samples (matched on sample_id and seq_type,
# or sample_id alone where the selection has no seq_type)
.selection_for = function(sel, sample_id, seq_type){
  i = match(paste(sample_id, seq_type), paste(sel$sample_id, sel$seq_type))
  untyped = which(is.na(sel$seq_type))
  j = untyped[match(sample_id, sel$sample_id[untyped])]
  sel[ifelse(is.na(i), j, i), , drop = FALSE]
}

# calibrate = TRUE in get_cn_segments: segments of samples from a relative
# copy number pipeline converted to absolute copy number with
# GAMBLR.utils::calibrate_relative_cn, against the same sample's
# segmentation from the first reference pipeline that passes QC (eligible,
# no contradiction penalty). Fits use genome-wide segments read from the
# per-pipeline files, and are applied to the rows returned (e.g. a region).
# Baselines of calibrated samples are recomputed from the calibrated CN.
.calibrate_relative_segments = function(segs, references, projection, baseline_method, max_CN, verbose){
  segs = as.data.frame(segs)
  if(!"pipeline" %in% colnames(segs)){
    message("calibrate: the segments have no pipeline column; nothing calibrated")
    return(list(segs = segs, calibration = NULL))
  }
  types = GAMBLR.utils::cn_pipeline_types()
  segs$cn_type = types$cn_type[match(segs$pipeline, types$pipeline)]
  # the relative CN, kept for relative pipelines' rows (calibrated or not)
  segs$CN_relative = ifelse(segs$cn_type %in% "relative", segs$CN, NA_real_)
  if(!"seg_seq_type" %in% colnames(segs)) segs$seg_seq_type = NA_character_
  relative = segs$cn_type %in% "relative" & !segs$pipeline %in% .published_cn_sources()
  todo = unique(segs[relative, c("ID", "pipeline", "seg_seq_type"), drop = FALSE])
  if(nrow(todo) == 0) return(list(segs = segs, calibration = NULL))
  no_baseline_fn = function(x, path, table = NULL) x
  read_segments = function(pipeline, seq_type, ids){
    s = tryCatch(suppressMessages(.get_pipeline_segments(pipeline, FALSE, projection, stats::setNames(list(ids), seq_type),
                                                         function(x) x, no_baseline_fn, verbose,
                                                         quiet_missing = TRUE, summarise = FALSE)),
                 error = function(e) NULL)
    if(is.null(s)) return(NULL)
    s = as.data.frame(s)
    if(!"log.ratio" %in% colnames(s) && "seg.mean" %in% colnames(s)) s$log.ratio = s$seg.mean
    if("dummy_segment" %in% colnames(s)) s = s[!s$dummy_segment %in% 1, , drop = FALSE]
    s
  }
  calibration = list()
  for(st in unique(todo$seg_seq_type)){
    in_st = todo[todo$seg_seq_type %in% st, , drop = FALSE]
    seq_type = if(is.na(st)) "capture" else st
    # each sample's reference: the first in `references` that passes QC
    scores = tryCatch(.cn_qc_scores(seq_type, projection, verbose), error = function(e) NULL)
    passing = if(is.null(scores)) NULL else
      scores[scores$pipeline %in% references & scores$eligible %in% TRUE & scores$pen_contradiction == 0, , drop = FALSE]
    in_st$reference = vapply(in_st$ID, function(id){
      ok = passing$pipeline[passing$sample_id == id]
      hit = references[references %in% ok]
      if(length(hit) == 0) NA_character_ else hit[1]
    }, "")
    ref_segs = do.call(rbind, lapply(stats::na.omit(unique(in_st$reference)), function(p){
      read_segments(p, seq_type, in_st$ID[in_st$reference %in% p])[, c("ID", "chrom", "start", "end", "CN")]
    }))
    for(p in unique(in_st$pipeline)){
      ids = in_st$ID[in_st$pipeline == p]
      full = read_segments(p, seq_type, ids)
      fits = if(is.null(full) || is.null(ref_segs)) NULL else
        attr(suppressMessages(GAMBLR.utils::calibrate_relative_cn(full, ref_segs)), "calibration")
      if(is.null(fits)) fits = data.frame(sample_id = ids, calibrated = FALSE, flag = "no_reference",
                                          stringsAsFactors = FALSE)
      fits = fits[match(ids, fits$sample_id), , drop = FALSE]
      fits$sample_id = ids
      fits$calibrated[is.na(fits$calibrated)] = FALSE
      fits$reference = in_st$reference[match(ids, in_st$ID)]
      fits$flag[is.na(fits$reference)] = "no_passing_reference"
      fits$calibrated[is.na(fits$reference)] = FALSE
      fits$pipeline = p
      fits$seq_type = st
      # apply each fit to the rows returned
      ok = fits[fits$calibrated, , drop = FALSE]
      rows = which(segs$pipeline == p & segs$seg_seq_type %in% st & segs$ID %in% ok$sample_id)
      at = match(segs$ID[rows], ok$sample_id)
      cn = pmax(0, (2^segs$log.ratio[rows] - ok$intercept[at]) / ok$slope[at])
      if(max_CN > 0) cn = pmin(cn, max_CN)
      segs$CN[rows] = round(cn, 2)
      segs$cn_type[rows] = "absolute"
      # baselines from the calibrated genome-wide CN
      if(length(rows) > 0 && "baseline_cn" %in% colnames(segs) && !is.null(baseline_method) && !is.null(full)){
        full_ok = full[full$ID %in% ok$sample_id, , drop = FALSE]
        f_at = match(full_ok$ID, ok$sample_id)
        full_ok$CN = pmax(0, (2^full_ok$log.ratio - ok$intercept[f_at]) / ok$slope[f_at])
        b = .cn_baselines_from_segments(full_ok, projection)
        value = b[[paste0("baseline_", baseline_method)]]
        value[!.baseline_ok(b)] = NA
        segs$baseline_cn[rows] = value[match(segs$ID[rows], b$ID)]
      }
      calibration[[paste(st, p)]] = fits
    }
  }
  calibration = dplyr::bind_rows(calibration)
  calibration = calibration[, c("sample_id", "seq_type", "pipeline", "reference", "calibrated",
                                intersect(c("purity", "slope", "intercept", "r_squared", "n_states",
                                            "altered_fraction", "cn2_fraction"), colnames(calibration)),
                                "flag"), drop = FALSE]
  if(all(is.na(segs$seg_seq_type))) segs$seg_seq_type = NULL
  done = table(calibration$reference[calibration$calibrated])
  message("calibrated ", sum(calibration$calibrated), " of ", nrow(calibration), " relative-CN sample(s)",
          if(length(done) > 0) paste0(" (against ", paste(names(done), done, collapse = ", "), ")") else "",
          if(any(!calibration$calibrated)) paste0("; ", sum(!calibration$calibrated),
                                                  " kept relative (see attr(<result>, \"calibration\")$flag)") else "")
  list(segs = segs, calibration = calibration)
}
