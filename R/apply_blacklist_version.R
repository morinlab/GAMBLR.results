#' @title Apply a Blacklist Version to SSM Calls.
#'
#' @description Make SSM calls from the current GAMBL results match an
#' earlier version of the slms-3 blacklists, by dropping the calls that
#' version's lists removed.
#'
#' @details The GAMBL merged and per-sample MAFs are deblacklisted by the
#' pipeline. In the 2026-11 audit, positions shown to be real hotspots
#' (e.g. FOXC1 S82N, HLA-B Y108*, CD70 splice) were taken off the clean
#' blacklists and the affected MAFs were regenerated, so the current results
#' keep calls the earlier lists dropped. `blacklist_version = "2022"`
#' restores the earlier behaviour, for results consistent with analyses done
#' before the audit.
#'
#' The removed positions are kept in one file per build, under
#' `resources$blacklist$removed$template` in the GAMBLR config
#' (`{projection}_removed_from_clean_blacklist.tsv`, next to the clean
#' blacklists). The pipeline applied only the lists of each sample's native
#' build (metadata `genome_build`) and made the other projection by liftover,
#' so a call is dropped only when its sample's native-build list dropped that
#' position. The position is translated into the coordinates of `projection`
#' using the other-build column of the file. As in the pipeline, a position is
#' dropped only if every list holding it gave it a count of at least
#' `drop_threshold`; any allele at the position is dropped.
#'
#' Samples without a `genome_build` in the metadata are left unchanged, with
#' a message.
#'
#' @param maf_df A data frame of SSM calls in MAF format from GAMBL results
#' (e.g. from [get_ssm_by_samples]). Needs Chromosome, Start_Position and
#' Tumor_Sample_Barcode. If it has a maf_seq_type column, samples are matched
#' to the metadata on sample_id and seq_type.
#' @param these_samples_metadata Metadata with sample_id, genome_build and
#' (optionally) seq_type for the samples in `maf_df`. Default: all of
#' [get_gambl_metadata].
#' @param projection Genome build of the coordinates in `maf_df`, grch37
#' (default) or hg38.
#' @param blacklist_version "2026-11" (default; the current lists, so nothing
#' is dropped) or "2022" (the lists before the 2026-11 audit).
#' @param drop_threshold The pipeline's blacklist count threshold. Default 4.
#' @param verbose Report how many calls were dropped. Default TRUE.
#'
#' @return `maf_df` without the calls the requested blacklist version drops.
#' A maf_data object is returned as maf_data.
#'
#' @import dplyr glue GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#'   meta = get_gambl_metadata() %>% dplyr::filter(seq_type == "capture")
#'   current = get_coding_ssm(these_samples_metadata = meta,
#'                            this_seq_type = "capture")
#'   pre_audit = apply_blacklist_version(current, meta,
#'                                       blacklist_version = "2022")
#'   # the same, in one call
#'   pre_audit = get_coding_ssm(these_samples_metadata = meta,
#'                              this_seq_type = "capture",
#'                              blacklist_version = "2022")
#' }
apply_blacklist_version = function(maf_df,
                                   these_samples_metadata,
                                   projection = "grch37",
                                   blacklist_version = "2026-11",
                                   drop_threshold = 4,
                                   verbose = TRUE){
  blacklist_version = .check_blacklist_version(blacklist_version)
  if(blacklist_version == "2026-11" || is.null(maf_df) || nrow(maf_df) == 0){
    return(maf_df)
  }
  if(!projection %in% c("grch37", "hg38")){
    stop("projection must be either grch37 or hg38")
  }
  needed = c("Chromosome", "Start_Position", "Tumor_Sample_Barcode")
  if(!all(needed %in% colnames(maf_df))){
    stop("maf_df needs the columns ", paste(needed, collapse = ", "),
         " to apply blacklist_version")
  }
  if(missing(these_samples_metadata) || is.null(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata()
  }
  if(!all(c("sample_id", "genome_build") %in% colnames(these_samples_metadata))){
    stop("these_samples_metadata needs sample_id and genome_build to apply blacklist_version")
  }

  # native build per sample (grch37, hs37d5, hg19-* -> grch37; hg38, hg38-* -> hg38)
  use_seq_type = "maf_seq_type" %in% colnames(maf_df) &&
    "seq_type" %in% colnames(these_samples_metadata)
  meta = as.data.frame(these_samples_metadata)
  meta = meta[!is.na(meta$genome_build) & meta$genome_build != "", , drop = FALSE]
  meta_key = if(use_seq_type) paste(meta$sample_id, meta$seq_type) else meta$sample_id
  meta_native = ifelse(grepl("38", meta$genome_build), "hg38", "grch37")
  if(!use_seq_type){
    ambiguous = unique(meta_key[duplicated(meta_key)])
    ambiguous = ambiguous[vapply(ambiguous, function(k) length(unique(meta_native[meta_key == k])) > 1, TRUE)]
    if(length(ambiguous) > 0){
      stop(length(ambiguous), " sample_id(s) have more than one native genome build in ",
           "these_samples_metadata (e.g. ", ambiguous[1], "); give metadata for one seq_type, ",
           "or a maf_df with a maf_seq_type column")
    }
  }
  row_key = if(use_seq_type){
    paste(maf_df$Tumor_Sample_Barcode, maf_df$maf_seq_type)
  }else{
    maf_df$Tumor_Sample_Barcode
  }
  native = meta_native[match(row_key, meta_key)]
  if(any(is.na(native))){
    unknown = unique(maf_df$Tumor_Sample_Barcode[is.na(native)])
    message(length(unknown), " sample(s) have no genome_build in these_samples_metadata ",
            "(e.g. ", unknown[1], "); their calls were kept unchanged by blacklist_version")
  }

  removed = .removed_blacklist_positions(projection, drop_threshold)
  chrom = sub("^chr", "", as.character(maf_df$Chromosome))
  drop = !is.na(native) &
    paste(native, chrom, format(as.numeric(maf_df$Start_Position), scientific = FALSE, trim = TRUE)) %in%
    paste(removed$native_build, removed$chrom, removed$pos)

  if(verbose){
    message("blacklist_version \"", blacklist_version, "\": dropped ", sum(drop), " call(s) from ",
            length(unique(maf_df$Tumor_Sample_Barcode[drop])), " sample(s) at positions removed ",
            "from the blacklists in the 2026-11 audit")
  }
  was_maf_data = inherits(maf_df, "maf_data")
  maf_df = maf_df[!drop, , drop = FALSE]
  if(was_maf_data){
    maf_df = GAMBLR.utils::create_maf_data(as.data.frame(maf_df), projection)
  }
  maf_df
}

.check_blacklist_version = function(blacklist_version){
  versions = c("2026-11", "2022")
  if(length(blacklist_version) != 1 || !blacklist_version %in% versions){
    stop("blacklist_version must be one of ", paste0('"', versions, '"', collapse = " or "))
  }
  blacklist_version
}

# Positions the pre-audit lists dropped, per native build, in the coordinates
# of `projection`: data frame of native_build, chrom (no "chr"), pos.
# A position counts only if every list row for it had count >= drop_threshold
# (a count below it on any list let the call through, as in the pipeline);
# an NA count (count not recorded) is taken as dropped.
.removed_blacklist_positions = function(projection, drop_threshold = 4){
  template = check_config_and_value("resources$blacklist$removed$template")
  base_path = check_config_and_value("project_base")
  builds = c("grch37", "hg38")
  out = list()
  for(native_build in builds){
    path = paste0(base_path, glue::glue(template, projection = native_build))
    if(!file.exists(path)){
      stop("cannot find the removed-blacklist positions file: ", path,
           "\nSync it from level_3 of the GAMBL results, or check resources$blacklist$removed$template in your config")
    }
    removed = utils::read.delim(path, colClasses = "character", check.names = FALSE)
    other_col = paste0(setdiff(builds, native_build), "_chrpos")
    if(!all(c("chrpos", "blacklist_count", other_col) %in% colnames(removed))){
      stop(path, " needs the columns chrpos, blacklist_count and ", other_col)
    }
    count = suppressWarnings(as.numeric(removed$blacklist_count))
    held = is.na(count) | count >= drop_threshold
    dropped = unique(removed$chrpos[removed$chrpos %in% removed$chrpos[held] &
                                      !removed$chrpos %in% removed$chrpos[!held]])
    if(projection == native_build){
      chrpos = dropped
    }else{
      chrpos = removed[[other_col]][match(dropped, removed$chrpos)]
      chrpos = chrpos[!is.na(chrpos) & chrpos != ""]
    }
    chrpos = sub("^chr", "", chrpos)
    out[[native_build]] = data.frame(native_build = native_build,
                                     chrom = sub(":.*$", "", chrpos),
                                     pos = sub("^.*:", "", chrpos))
  }
  do.call(rbind, out)
}
