#' @title Get Battenberg output files
#'
#' @description Find the Battenberg output files (allele-specific copy
#' number, purity/ploidy and fit details) for tumour genomes with a matched
#' normal, across unix groups (e.g. gambl and icgc_dart), using the templates
#' in the GAMBLR config. Returns a table of file paths, one row per sample,
#' for use with [GAMBLR.utils::summarize_battenberg] or for reading the files
#' directly.
#'
#' @details Paths are `project_base` (from the config, or the
#' `project_base` argument) plus the `results_flatfiles$cnv` templates
#' `battenberg_dir` and `battenberg_files` (built-in defaults are used if the
#' config doesn't define them), filled in from each
#' sample's metadata (`unix_group`, `seq_type`, `genome_build`, `sample_id`,
#' `normal_sample_id`). Battenberg is only run on genomes with a matched
#' normal (`battenberg_wildcards`), so other samples are left out with a
#' message. `projection` is the coordinate system of the files: grch37 for
#' grch37/hs37d5/hg19 alignments, hg38 for hg38/grch38.
#'
#' @param these_samples_metadata Metadata with sample_id, normal_sample_id,
#' unix_group, seq_type, genome_build and pairing_status, e.g. from
#' [GAMBLR.results::get_gambl_metadata]. Default: all genomes from
#' `get_gambl_metadata()`.
#' @param file_types Which files to locate (names in the config's
#' `battenberg_files`). Default: subclones, cellularity_ploidy, rho_and_psi
#' and solution_status.
#' @param project_base Optional base directory replacing the config's
#' `project_base`, e.g. a local copy of the results that keeps the same
#' directory layout (`<project_base>/<unix_group>/battenberg-1.1/...`).
#' @param require_all Set to TRUE to keep only samples for which every
#' requested file exists. Default FALSE (all eligible samples, with
#' `<type>_exists` columns).
#'
#' @return A data frame with sample_id, normal_sample_id, unix_group,
#' seq_type, genome_build, projection, battenberg_dir, and for each file
#' type a path column (`<type>_file`) and an `<type>_exists` column.
#'
#' @import dplyr glue GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' bb_files = get_battenberg(get_gambl_metadata(seq_type_filter = "genome"))
#' table(bb_files$unix_group, bb_files$subclones_exists)
#'
#' # summarize all available samples
#' bb = GAMBLR.utils::summarize_battenberg(bb_files[bb_files$subclones_exists, ])
#' }
get_battenberg = function(these_samples_metadata,
                          file_types = c("subclones", "cellularity_ploidy",
                                         "rho_and_psi", "solution_status"),
                          project_base = NULL,
                          require_all = FALSE){
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome")
  }
  required = c("sample_id", "normal_sample_id", "unix_group", "seq_type",
               "genome_build", "pairing_status")
  missing_cols = setdiff(required, colnames(these_samples_metadata))
  if(length(missing_cols) > 0){
    stop(paste("these_samples_metadata is missing columns:",
               paste(missing_cols, collapse = ", "),
               "(use GAMBLR.results::get_gambl_metadata())"))
  }

  base = if(is.null(project_base)){
    GAMBLR.helpers::check_config_value(config::get("project_base"))
  }else{
    sub("/*$", "/", project_base)
  }
  # templates from the config, or these defaults if the config (e.g. an older
  # config.yml in the working directory) doesn't define them
  cnv_config = tryCatch(config::get("results_flatfiles")$cnv,
                        error = function(e) NULL)
  dir_template = cnv_config$battenberg_dir
  file_templates = cnv_config$battenberg_files
  if(is.null(dir_template) || is.null(file_templates)){
    message("battenberg_dir/battenberg_files not in the config; using the default templates")
    dir_template = "{unix_group}/battenberg-1.1/02-battenberg/{seq_type}--{genome_build}/{sample_id}--{normal_sample_id}/"
    file_templates = list(subclones = "{sample_id}_subclones.txt",
                          cellularity_ploidy = "{sample_id}_cellularity_ploidy.txt",
                          rho_and_psi = "{sample_id}_rho_and_psi.txt",
                          solution_status = "{sample_id}_copynumber_solution_status.txt")
  }
  unknown = setdiff(file_types, names(file_templates))
  if(length(unknown) > 0){
    stop(paste0("unknown file_types: ", paste(unknown, collapse = ", "),
                "; the config defines: ",
                paste(names(file_templates), collapse = ", ")))
  }
  pairing = if(!is.null(cnv_config$battenberg_wildcards))
    cnv_config$battenberg_wildcards else list(pairing_status = "matched")

  # Battenberg runs on genomes with a matched normal
  meta = as.data.frame(these_samples_metadata)
  eligible = meta$seq_type == "genome" &
    meta$pairing_status %in% unname(unlist(pairing)) &
    !is.na(meta$normal_sample_id) & meta$normal_sample_id != ""
  if(any(!eligible)){
    message(sum(!eligible), " sample(s) are not genomes with a matched normal ",
            "and have no Battenberg output; they are left out")
  }
  meta = meta[eligible, , drop = FALSE]
  if(nrow(meta) == 0){
    stop("no genomes with a matched normal in these_samples_metadata")
  }

  out = data.frame(sample_id = meta$sample_id,
                   normal_sample_id = meta$normal_sample_id,
                   unix_group = meta$unix_group,
                   seq_type = meta$seq_type,
                   genome_build = meta$genome_build,
                   stringsAsFactors = FALSE)
  build = tolower(out$genome_build)
  out$projection = ifelse(build %in% c("grch37", "hs37d5", "hg19"), "grch37",
                          ifelse(build %in% c("hg38", "grch38"), "hg38",
                                 NA_character_))
  if(anyNA(out$projection)){
    warning("unrecognised genome_build for ", sum(is.na(out$projection)),
            " sample(s): ", paste(unique(out$genome_build[is.na(out$projection)]),
                                  collapse = ", "))
  }
  fill = function(template){
    as.character(glue::glue_data(out, template))
  }
  out$battenberg_dir = paste0(base, fill(dir_template))
  for(type in file_types){
    # battenberg_dir ends with "/"
    path = paste0(sub("/*$", "/", out$battenberg_dir), fill(file_templates[[type]]))
    out[[paste0(type, "_file")]] = path
    out[[paste0(type, "_exists")]] = file.exists(path)
  }

  # report what was found, per unix group
  found = out[[paste0(file_types[1], "_exists")]]
  counts = table(out$unix_group, factor(found, levels = c(TRUE, FALSE)))
  message("Battenberg ", file_types[1], " files found: ",
          paste0(rownames(counts), " ", counts[, "TRUE"], "/",
                 rowSums(counts), collapse = "; "))
  if(require_all){
    all_found = Reduce(`&`, lapply(file_types, function(type)
      out[[paste0(type, "_exists")]]))
    out = out[all_found, , drop = FALSE]
  }
  rownames(out) = NULL
  out
}
