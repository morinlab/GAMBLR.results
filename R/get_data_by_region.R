#' @title Get mutations, SVs and copy number segments for a region
#'
#' @description Retrieve the data used by GAMBLR.viz region plots
#' (`plot_ssm_density`, `plot_sv_density`, `plot_cohort_CN_region`) for one
#' region in one call: simple somatic mutations in the region, structural
#' variants with a breakpoint in it, and copy number segments, for the
#' samples in `these_samples_metadata`, together with that metadata
#' annotated with which data types are available for each sample.
#'
#' @details Calls [get_ssm_by_region] (full MAF columns), [get_combined_sv]
#' (or [get_manta_sv]) and [get_cn_segments]. SVs are kept when either
#' breakpoint is in the region, so translocations into it are included.
#'
#' \strong{Availability.} A region query can't tell a sample with no events
#' in the region from a sample with no data, so the returned metadata (all
#' rows kept) gets `ssm_available`, `sv_available`, `cn_available` and
#' `all_available` columns, as from [annotate_data_availability]. Each data
#' type is returned only for samples where it is available. GAMBLR.viz
#' region plots use these columns, when present, for their denominators and
#' exclusions.
#'
#' A sample_id with more than one copy number profile (genome and capture)
#' keeps only the profile for the seq_type of its row in
#' `these_samples_metadata`. Copy number segments are those overlapping the
#' region, each with its sample's genome-wide baseline copy number
#' (`baseline_cn`, `baseline_method`), from a cache built once per seg file.
#'
#' @param region Region as "chrom:start-end", e.g. "chr3:187,425,881-187,500,524".
#' @param these_samples_metadata Metadata (e.g. from [get_gambl_metadata]);
#' all data are restricted to its samples.
#' @param projection "grch37" (default) or "hg38".
#' @param data_types Which data to retrieve: any of "ssm", "sv" and "cn".
#' Default all three.
#' @param sv_source "combined" (default; [get_combined_sv]) or "manta"
#' ([get_manta_sv]).
#' @param baseline_method Baseline copy number added to the segments: "mode"
#' (default), "mode_rounded" or "mean" (see
#' [GAMBLR.utils::add_segment_baseline]).
#' @param ... Other arguments passed to [get_ssm_by_region] (e.g.
#' `min_read_support`).
#'
#' @return A list with `region`, `metadata` (`these_samples_metadata` with
#' the availability columns), and for the requested data types `ssm` (MAF
#' rows in the region), `sv` (SVs with a breakpoint in the region), `seg`
#' (copy number segments in the region, with baselines), and `ssm_samples`
#' and `sv_samples` (IDs of the samples with that data type).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(pathology == "DLBCL")
#' d = get_data_by_region("chr3:187,425,881-187,500,524", meta)
#' table(d$metadata$all_available)
#'
#' GAMBLR.viz::plot_cohort_CN_region(d$seg, region = d$region, genes = "BCL6",
#'                                   these_samples_metadata = d$metadata,
#'                                   split_by = "lymphgen",
#'                                   sv_data = d$sv, ssm_data = d$ssm)
#' }
get_data_by_region = function(region,
                              these_samples_metadata,
                              projection = "grch37",
                              data_types = c("ssm", "sv", "cn"),
                              sv_source = "combined",
                              baseline_method = "mode",
                              ...){
  data_types = match.arg(data_types, c("ssm", "sv", "cn"), several.ok = TRUE)
  sv_source = match.arg(sv_source, c("combined", "manta"))
  if(missing(these_samples_metadata)){
    stop("these_samples_metadata is required; the data are restricted to its samples")
  }
  # parse the region once (commas allowed, chr prefix optional)
  clean = gsub(",", "", region)
  parts = regmatches(clean, regexec("^(chr)?([0-9XY]+):([0-9.e+]+)-([0-9.e+]+)$", clean))[[1]]
  if(length(parts) == 0){
    stop("region must look like \"chr3:187425881-187500524\"")
  }
  chrom = parts[3]
  start = as.numeric(parts[4]); end = as.numeric(parts[5])
  region = paste0(chrom, ":", format(start, scientific = FALSE), "-",
                  format(end, scientific = FALSE))
  meta = these_samples_metadata
  out = list(region = region)

  if("ssm" %in% data_types){
    meta = annotate_ssm_availability(meta, projection = projection)
    ssm_meta = meta[meta$ssm_available, , drop = FALSE]
    ssm = get_ssm_by_region(region = region, these_samples_metadata = ssm_meta,
                            projection = projection, streamlined = FALSE, ...)
    out$ssm = ssm[ssm$Tumor_Sample_Barcode %in% ssm_meta$sample_id, , drop = FALSE]
    out$ssm_samples = unique(ssm_meta$sample_id)
  }
  if("sv" %in% data_types){
    meta = annotate_sv_availability(meta, projection = projection, sv_source = sv_source)
    sv_meta = meta[meta$sv_available, , drop = FALSE]
    if(nrow(sv_meta) > 0){
      sv = if(sv_source == "combined"){
        get_combined_sv(these_samples_metadata = sv_meta, projection = projection, region = region)
      }else{
        get_manta_sv(these_samples_metadata = sv_meta, projection = projection, region = region)
      }
      out$sv = sv[sv$tumour_sample_id %in% sv_meta$sample_id, , drop = FALSE]
    }
    out$sv_samples = unique(sv_meta$sample_id)
  }
  if("cn" %in% data_types){
    meta = annotate_cn_availability(meta, projection = projection)
    cn_meta = meta[meta$cn_available, , drop = FALSE]
    # region segments, each with its sample's genome-wide baseline
    seg = get_cn_segments(these_samples_metadata = cn_meta, projection = projection,
                          region = region, baseline_method = baseline_method)
    # one profile per sample: the seq_type of its metadata row (a sample_id
    # can have both genome and capture profiles, which would double count)
    if("seg_seq_type" %in% colnames(seg)){
      wanted = paste(cn_meta$sample_id, cn_meta$seq_type)
      keep = is.na(seg$seg_seq_type) | paste(seg$ID, seg$seg_seq_type) %in% wanted
      seg = seg[keep, , drop = FALSE]
    }
    out$seg = seg
  }
  out$metadata = meta
  out
}
