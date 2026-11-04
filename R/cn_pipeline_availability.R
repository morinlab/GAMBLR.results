#' @title Which copy number pipelines have segments for each sample
#'
#' @description For the samples in `these_samples_metadata`, count each
#' pipeline's segments in the per-pipeline seg files (as used by
#' [get_cn_segments] with `pipeline`), and report the module the combined
#' merge uses for each sample. Use it to see what can be compared for a
#' sample before loading the segments, e.g. with
#' `get_cn_segments(..., pipeline = "all")`.
#'
#' @details Only the ID column of each per-pipeline file (and ID and module of
#' the combined file) is read, so this is quick even for genome merges.
#' Dummy (filled) segments are counted where the files still contain them;
#' [get_cn_segments] drops them.
#'
#' @param these_samples_metadata Metadata with sample_id and seq_type (genome
#' and capture samples are used).
#' @param projection Genome build ("grch37" or "hg38"). Default "grch37".
#' @param pipeline Pipelines to check. Default "all": every pipeline with a
#' per-pipeline seg file for the seq_types in the metadata, and the published
#' sources in GAMBLR.data (for capture samples).
#'
#' @return A data frame with one row per sample: sample_id, seq_type, one
#' column per pipeline with its number of segments for the sample (0 when its
#' file lacks the sample, NA when the pipeline has no file for the sample's
#' seq_type), n_pipelines (pipelines with segments), and combined_module (the
#' module of the sample's segments in the combined merge, NA if it isn't
#' there).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(cohort == "dlbcl_chapuy")
#' cn_pipeline_availability(meta)
#' }
cn_pipeline_availability = function(these_samples_metadata, projection = "grch37",
                                    pipeline = "all"){
  meta = as.data.frame(these_samples_metadata)
  meta = meta[meta$seq_type %in% c("genome", "capture"), c("sample_id", "seq_type"), drop = FALSE]
  if(nrow(meta) == 0) stop("no genome or capture samples in these_samples_metadata")
  seq_types = unique(meta$seq_type)
  pipelines = .resolve_cn_pipelines(pipeline, NULL, seq_types, projection)$pipeline

  read_ids = function(path, cols = "ID"){
    suppressMessages(readr::read_tsv(path, col_select = dplyr::all_of(cols),
                                     col_types = readr::cols(.default = "c"),
                                     progress = FALSE, lazy = TRUE))
  }
  out = do.call(rbind, lapply(seq_types, function(st){
    m = meta[meta$seq_type == st, , drop = FALSE]
    for(p in pipelines){
      m[[p]] = if(.pipeline_has_data(p, st, projection)){
        ids = if(p %in% .published_cn_sources()) .published_segments(p, projection)$ID else
          read_ids(.cnv_pipeline_path(p, st, projection))$ID
        counts = table(ids)
        n = as.integer(counts[m$sample_id])
        ifelse(is.na(n), 0L, n)
      } else NA_integer_
    }
    combined = .cnv_combined_path(st, projection)
    m$combined_module = if(file.exists(combined)){
      header = colnames(suppressMessages(readr::read_tsv(combined, n_max = 0, progress = FALSE)))
      if("module" %in% header){
        cm = read_ids(combined, c("ID", "module"))
        cm = cm[!duplicated(cm$ID), , drop = FALSE]
        cm$module[match(m$sample_id, cm$ID)]
      } else NA_character_
    } else NA_character_
    m
  }))
  out$n_pipelines = rowSums(as.matrix(out[, pipelines, drop = FALSE]) > 0, na.rm = TRUE)
  out = out[, c("sample_id", "seq_type", pipelines, "n_pipelines", "combined_module")]
  rownames(out) = NULL
  out
}
