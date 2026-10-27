#' @title Annotate metadata with the data types available for each sample
#'
#' @description Add `*_available` columns to a metadata table saying, for
#' each sample, whether it has simple somatic mutation (`ssm_available`),
#' structural variant (`sv_available`) and copy number (`cn_available`) data,
#' plus `all_available` (whichever of these three columns are present are
#' all TRUE; other columns such as `bam_available` are not considered).
#' All rows are kept unless `drop_unavailable` is TRUE, so the table still
#' shows the full cohort and which samples each data type covers.
#'
#' @details A query of a merged results file can't tell a sample with no
#' events (in a region, gene or otherwise) from a sample that isn't in the
#' file at all, so availability is determined separately:
#' \itemize{
#'   \item `ssm_available`: the sample is in the MAF for its seq_type (the
#'     merged, tabix-indexed MAF read by [get_ssm_by_regions], or with
#'     `source = "cds"` the coding MAF read by [get_coding_ssm]), from a
#'     sample index built once per file (see [get_file_samples]). For the
#'     merged MAF, samples listed by [get_excluded_samples] are FALSE, since
#'     the region getters skip them. Capture samples only cover their
#'     targets, so TRUE doesn't mean a capture sample covers a non-coding
#'     region.
#'   \item `sv_available`: a genome present in the merged SV file of
#'     `sv_source` ("combined" as read by [get_combined_sv], or "manta" as
#'     read by [get_manta_sv]).
#'   \item `cn_available`: the sample has copy number segments (combined seg
#'     file) that cover enough of its genome for a baseline and have no
#'     coordinates past the chromosome ends; see [get_cn_baselines].
#' }
#' Calling one of these on a table that already has other `*_available`
#' columns keeps them and updates `all_available`, so they can be chained.
#' The getters take `with_metadata = TRUE` to return their data together
#' with the metadata annotated this way.
#'
#' Building a sample index or baseline cache for a file the first time
#' takes a while (minutes for the genome MAF); after that it is read from
#' the cache until the results file changes.
#'
#' @param these_samples_metadata Metadata with sample_id and seq_type.
#' @param data_types Any of "ssm", "sv" and "cn". Default all three.
#' @param projection "grch37" (default) or "hg38".
#' @param source For SSMs, the MAF the data come from: "merged" (default;
#' the indexed MAF used by the region getters) or "cds" (the coding MAF).
#' @param augmented For SSMs, the augmented (default) or deblacklisted MAF.
#' @param seq_types For SSMs, the seq_types the data come from; rows of
#' other seq_types are FALSE. Default genome and capture.
#' @param sv_source For SVs, "combined" (default) or "manta".
#' @param pipeline,priority For copy number, optionally the per-pipeline seg
#' files to use instead of the combined one, as in [get_cn_segments]. The
#' metadata then also gets `cn_pipeline`, the pipeline each sample's copy
#' number would come from.
#' @param min_genome_coverage,min_autosomes For copy number, as in
#' [get_cn_baselines]. Defaults 0.5 and 20.
#' @param drop_unavailable Set to TRUE to also drop the rows that lack the
#' data type (for `annotate_data_availability`, rows with `all_available`
#' FALSE). Default FALSE: rows are only annotated.
#' @param verbose Report how many samples lack the data type. Default TRUE.
#'
#' @return `these_samples_metadata` with the `*_available` columns and
#' `all_available`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(pathology == "DLBCL")
#' meta = annotate_data_availability(meta)
#' table(meta$seq_type, meta$all_available)
#'
#' # one data type, chained with another
#' meta = annotate_ssm_availability(meta) %>% annotate_cn_availability()
#' }
annotate_data_availability = function(these_samples_metadata,
                                      data_types = c("ssm", "sv", "cn"),
                                      projection = "grch37",
                                      sv_source = "combined",
                                      drop_unavailable = FALSE,
                                      verbose = TRUE){
  data_types = match.arg(data_types, c("ssm", "sv", "cn"), several.ok = TRUE)
  meta = these_samples_metadata
  if("ssm" %in% data_types){
    meta = annotate_ssm_availability(meta, projection = projection, verbose = verbose)
  }
  if("sv" %in% data_types){
    meta = annotate_sv_availability(meta, projection = projection, sv_source = sv_source,
                                    verbose = verbose)
  }
  if("cn" %in% data_types){
    meta = annotate_cn_availability(meta, projection = projection, verbose = verbose)
  }
  if(drop_unavailable){
    meta = meta[meta$all_available, , drop = FALSE]
  }
  meta
}

#' @rdname annotate_data_availability
#' @export
annotate_ssm_availability = function(these_samples_metadata,
                                     projection = "grch37",
                                     source = "merged",
                                     augmented = TRUE,
                                     seq_types = c("genome", "capture"),
                                     drop_unavailable = FALSE,
                                     verbose = TRUE){
  source = match.arg(source, c("merged", "cds"))
  meta = these_samples_metadata
  available = rep(FALSE, nrow(meta))
  for(st in intersect(unique(meta$seq_type), intersect(seq_types, c("genome", "capture")))){
    rows = meta$seq_type %in% st
    ids = .file_samples_or_none(.ssm_file_path(st, projection, source, augmented),
                                "Tumor_Sample_Barcode")
    available[rows] = meta$sample_id[rows] %in% ids
  }
  if(source == "merged"){
    # the region getters skip these
    available[meta$sample_id %in% get_excluded_samples("slms_3")] = FALSE
  }
  .set_available(meta, "ssm_available", available, "SSM data", drop_unavailable, verbose)
}

#' @rdname annotate_data_availability
#' @export
annotate_sv_availability = function(these_samples_metadata,
                                    projection = "grch37",
                                    sv_source = "combined",
                                    drop_unavailable = FALSE,
                                    verbose = TRUE){
  sv_source = match.arg(sv_source, c("combined", "manta"))
  meta = these_samples_metadata
  ids = .file_samples_or_none(.sv_merged_path(sv_source, projection), "tumour_sample_id")
  available = meta$seq_type %in% "genome" & meta$sample_id %in% ids
  .set_available(meta, "sv_available", available, "SV data", drop_unavailable, verbose)
}

#' @rdname annotate_data_availability
#' @export
annotate_cn_availability = function(these_samples_metadata,
                                    projection = "grch37",
                                    pipeline = NULL,
                                    priority = NULL,
                                    min_genome_coverage = 0.5,
                                    min_autosomes = 20,
                                    drop_unavailable = FALSE,
                                    verbose = TRUE){
  meta = these_samples_metadata
  b = get_cn_baselines(meta, projection = projection, pipeline = pipeline, priority = priority,
                       min_genome_coverage = min_genome_coverage, min_autosomes = min_autosomes)
  at = match(paste(meta$sample_id, meta$seq_type), paste(b$sample_id, b$seq_type))
  available = b$ok[at] %in% TRUE
  if(!is.null(pipeline)) meta$cn_pipeline = b$pipeline[at]
  .set_available(meta, "cn_available", available, "copy number data", drop_unavailable, verbose)
}

# Set one availability column, update all_available, report and optionally
# drop the rows without the data type
.set_available = function(meta, column, available, what, drop_unavailable, verbose){
  meta[[column]] = available
  columns = intersect(c("ssm_available", "sv_available", "cn_available"), colnames(meta))
  meta$all_available = Reduce(`&`, lapply(columns, function(col) meta[[col]] %in% TRUE))
  if(verbose && any(!available)){
    message(sum(!available), " of ", nrow(meta), " samples have no usable ", what,
            " (", column, " is FALSE)")
  }
  if(drop_unavailable) meta = meta[available, , drop = FALSE]
  meta
}

# Sample IDs in a results file (from its cached index), or none with a
# warning when the file isn't there
.file_samples_or_none = function(file, sample_column){
  if(!file.exists(file)){
    warning("results file not found, so no samples are marked available from it: ", file,
            call. = FALSE)
    return(character(0))
  }
  get_file_samples(file, sample_column = sample_column)
}

# Path of the MAF read by the SSM getters for a seq_type: the merged
# (tabix-indexed) MAF or the coding (cds) MAF
.ssm_file_path = function(seq_type, projection, source = "merged", augmented = TRUE){
  key = paste0("results_flatfiles$ssm$template$", source, "$",
               if(augmented) "augmented" else "deblacklisted")
  path = paste0(check_config_and_value("project_base"), glue::glue(check_config_and_value(key)))
  if(source == "merged") path = paste0(path, ".bgz")
  path
}

# Path of the merged SV file used by get_combined_sv() or get_manta_sv()
.sv_merged_path = function(sv_source, projection){
  base = check_config_and_value("project_base")
  if(sv_source == "combined"){
    f = check_config_and_value("results_flatfiles$sv_combined$icgc_dart")
    if(projection == "hg38") f = sub("--grch37", "--hg38", f)
    path = paste0(base, f)
    if(file.access(path, 4) == -1){
      path = paste0(base, check_config_and_value("results_flatfiles$sv_combined$gambl"))
    }
    return(path)
  }
  path = paste0(base, glue::glue(check_config_and_value("results_merged$manta_sv$icgc_dart")))
  if(file.access(path, 4) == -1){
    path = paste0(base, glue::glue(check_config_and_value("results_merged$manta_sv$gambl")))
  }
  path
}

# For getters called with with_metadata = TRUE: run the same call without
# it and return the data with the metadata annotated by `annotate`
.with_metadata = function(call, env, data_name, annotate){
  call$with_metadata = FALSE
  out = list(eval(call, env))
  names(out) = data_name
  out$metadata = annotate()
  out
}
