#' @title Compute per-sample mutation burden for a sample scope.
#'
#' @description Core computation for the `mutation_burden` table of
#' [collate_results_db()]: counts each sample's SLMS-3 calls (pooled, SNVs
#' and indels) and a burden per Mb over a territory that depends on the
#' assay:
#' \itemize{
#'   \item genomes: every call, over the effective genome size of the build
#'     (`mutation_burden$effective_genome_mb` in the config);
#'   \item exomes with a known kit: calls inside the kit's unpadded target
#'     regions (`mutation_burden$capture_space_bed`, merged), over the
#'     targets' size;
#'   \item exomes whose kit isn't known (`mutation_burden$fallback_capture_spaces`,
#'     or a kit whose target BED can't be found): calls in coding sequence
#'     (`burden_cds_classes`, which include Silent), over the merged
#'     protein-coding CDS of the build (`mutation_burden$fallback_coding_exome_mb`).
#' }
#' Counts are read in the projection of each sample's own build (grch37 or
#' hg38), so kit targets and calls share coordinates. A call is inside a
#' target when its start position is.
#'
#' @param sample_table A data frame with sample_id, seq_type, genome_build,
#' capture_space and the columns `get_ssm_by_samples()` needs.
#'
#' @return A data frame with `sample_id`; `n_ssm`, `n_snv`, `n_indel` (all
#' calls); `burden_territory` ("genome", "capture_targets" or
#' "coding_fallback"), `burden_territory_source`, `burden_mb`;
#' `burden_mutations`, `burden_snv`, `burden_indel` (calls in the
#' territory); and `burden_per_mb`, `burden_snv_per_mb`,
#' `burden_indel_per_mb`.
#'
#' @import dplyr
#' @keywords internal
#' @noRd
compute_mutation_burden_core <- function(sample_table){
  settings = burden_settings()
  sample_table = as.data.frame(sample_table)
  sample_table$.family = burden_build_family(sample_table$genome_build)
  out = list()
  for(fam in unique(sample_table$.family)){
    rows = sample_table[sample_table$.family == fam, , drop = FALSE]
    # samples whose MAF doesn't exist get NA rows (from the join below) without being read
    readable = ssm_samples_with_input(rows, projection = fam)
    if(nrow(readable) == 0) next
    rows = readable
    muts = get_ssm_by_samples(
      these_samples_metadata = rows,
      projection = fam,
      augmented = FALSE,
      min_read_support = 0,
      basic_columns = TRUE
    )
    muts = as.data.frame(muts)
    muts = muts[muts$Tumor_Sample_Barcode %in% rows$sample_id, , drop = FALSE]
    muts = muts[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "Variant_Type",
                    "Variant_Classification")]
    for(i in seq_len(nrow(rows))){
      r = rows[i, , drop = FALSE]
      m = muts[muts$Tumor_Sample_Barcode == r$sample_id, , drop = FALSE]
      territory = burden_territory(r, fam, settings)
      in_territory = switch(territory$type,
        genome = rep(TRUE, nrow(m)),
        capture_targets = burden_in_targets(m, territory$targets),
        coding_fallback = m$Variant_Classification %in% burden_cds_classes)
      snv = m$Variant_Type %in% "SNP"
      indel = m$Variant_Type %in% c("INS", "DEL")
      n_in = sum(in_territory)
      out[[length(out) + 1]] = data.frame(
        sample_id = r$sample_id,
        n_ssm = nrow(m), n_snv = sum(snv), n_indel = sum(indel),
        burden_territory = territory$type, burden_territory_source = territory$source,
        burden_mb = territory$mb,
        burden_mutations = n_in, burden_snv = sum(in_territory & snv), burden_indel = sum(in_territory & indel),
        burden_per_mb = round(n_in / territory$mb, 4),
        burden_snv_per_mb = round(sum(in_territory & snv) / territory$mb, 4),
        burden_indel_per_mb = round(sum(in_territory & indel) / territory$mb, 4),
        stringsAsFactors = FALSE)
    }
  }
  result = dplyr::bind_rows(out)
  if(nrow(result) == 0){
    result = data.frame(sample_id = character(0), n_ssm = integer(0), n_snv = integer(0), n_indel = integer(0),
                        burden_territory = character(0), burden_territory_source = character(0),
                        burden_mb = numeric(0), burden_mutations = integer(0), burden_snv = integer(0),
                        burden_indel = integer(0), burden_per_mb = numeric(0), burden_snv_per_mb = numeric(0),
                        burden_indel_per_mb = numeric(0), stringsAsFactors = FALSE)
  }
  # every requested sample gets a row, as for compute_ssm_results_core()
  dplyr::left_join(dplyr::select(as.data.frame(sample_table), sample_id), result, by = "sample_id")
}

# Coding-sequence variant classes counted for exomes of unknown kit (inside the CDS; splice and targeted-region
# classes can fall outside it)
burden_cds_classes = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                       "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Translation_Start_Site")

# grch37 or hg38: the projection matching a sample's build (e.g. hs37d5, hg19-reddy -> grch37; hg38-nci -> hg38)
burden_build_family = function(genome_build){
  ifelse(grepl("38", genome_build), "hg38", "grch37")
}

burden_settings = function(){
  s = check_config_and_value("mutation_burden")
  s$fallback_capture_spaces = trimws(strsplit(s$fallback_capture_spaces, ",")[[1]])
  s
}

# The territory of one sample: its type, size in Mb, a description of its source and, for kits, the merged targets
burden_territory = function(row, family, settings){
  if(row$seq_type == "genome"){
    return(list(type = "genome", mb = as.numeric(settings$effective_genome_mb[[family]]),
                source = paste0(settings$effective_genome_source, ", ", family)))
  }
  fallback = list(type = "coding_fallback", mb = as.numeric(settings$fallback_coding_exome_mb[[family]]),
                  source = paste0(settings$fallback_coding_exome_source, ", ", family))
  kit = row$capture_space
  if(is.na(kit) || !nzchar(kit) || kit %in% settings$fallback_capture_spaces) return(fallback)
  targets = burden_kit_targets(row$genome_build, kit, settings)
  if(is.null(targets)){
    fallback$source = paste0(fallback$source, " (no target BED for ", kit, " at ",
                             burden_kit_bed_path(row$genome_build, kit, settings), ")")
    return(fallback)
  }
  list(type = "capture_targets", mb = round(attr(targets, "bp") / 1e6, 4), targets = targets,
       source = paste0(kit, " unpadded targets (", attr(targets, "path"), ")"))
}

# A kit's unpadded target BED, merged, as one data frame per chromosome (chr prefix removed); NULL if missing.
# Read once per session per build and kit.
.burden_target_cache = new.env(parent = emptyenv())
burden_kit_targets = function(genome_build, capture_space, settings){
  key = paste(genome_build, capture_space, sep = "|")
  if(base::exists(key, envir = .burden_target_cache, inherits = FALSE)) return(base::get(key, envir = .burden_target_cache))
  path = burden_kit_bed_path(genome_build, capture_space, settings)
  targets = NULL
  if(file.exists(path)){
    bed = data.table::fread(path, header = FALSE, select = 1:3, col.names = c("chrom", "start", "end"))
    bed$chrom = sub("^chr", "", bed$chrom)
    bed = bed[order(bed$chrom, bed$start), ]
    merged = lapply(split(bed, bed$chrom), function(d){
      s = d$start; e = d$end
      new_run = c(TRUE, s[-1] > cummax(e)[-length(e)])
      run = cumsum(new_run)
      data.frame(start = tapply(s, run, min), end = tapply(e, run, max))
    })
    targets = merged
    attr(targets, "bp") = sum(vapply(merged, function(d) sum(d$end - d$start), numeric(1)))
    attr(targets, "path") = path
  }
  base::assign(key, targets, envir = .burden_target_cache)
  targets
}

burden_kit_bed_path = function(genome_build, capture_space, settings){
  paste0(check_config_and_value("repo_base"), glue::glue(settings$capture_space_bed))
}

# Which calls start inside the merged targets (BED: 0-based, half-open)
burden_in_targets = function(muts, targets){
  inside = logical(nrow(muts))
  chrom = sub("^chr", "", as.character(muts$Chromosome))
  for(ch in intersect(unique(chrom), names(targets))){
    t = targets[[ch]]
    sel = which(chrom == ch)
    pos0 = muts$Start_Position[sel] - 1
    idx = findInterval(pos0, t$start)
    inside[sel] = idx > 0 & pos0 < t$end[pmax(idx, 1)]
  }
  inside
}

#' @title Input MAF of each sample for the mutation burden, in its own build's projection.
#' @keywords internal
#' @noRd
mutation_burden_input_files = function(sample_table){
  sample_table = as.data.frame(sample_table)
  family = burden_build_family(sample_table$genome_build)
  out = sample_table[, "sample_id", drop = FALSE]
  out$input_path = NA_character_
  for(fam in unique(family)){
    sel = family == fam
    out$input_path[sel] = ssm_input_files(sample_table[sel, , drop = FALSE], projection = fam)$input_path
  }
  out
}
