#' @title Compute curated SV annotations for a sample scope.
#'
#' @description Core computation behind `collate_curated_sv_results()`,
#' extracted so it can also be called directly by [collate_results_db()]
#' for just the subset of samples missing from its cache table. Every
#' manually curated `.tsv` under the configured `derived_and_curated`
#' directory is left-joined on, preserving the original's implicit-`by`
#' behaviour (whichever columns each file happens to share with the
#' accumulating table -- typically `sample_id` or `biopsy_id`). Because the
#' join always starts from the full `sample_table` and uses `left_join`,
#' every requested `sample_id` keeps its row even if none of the curated
#' files had anything for it (NA in the new columns) -- only the columns
#' not already present on `sample_table` are kept in the return value, so
#' the cache table doesn't end up storing a duplicate copy of whatever
#' metadata was used purely as a join key.
#'
#' @param sample_table A data frame with sample_id as a column, scoping
#' which samples to compute for. Should include whatever metadata columns
#' (e.g. `biopsy_id`) the curated files are expected to join on.
#'
#' @return A data frame with `sample_id` plus whatever new columns the
#' curated files contribute.
#'
#' @import readr dplyr GAMBLR.helpers
#'
#' @keywords internal
#' @noRd
compute_curated_sv_results_core <- function(sample_table){

  path_to_files = GAMBLR.helpers::check_config_value(config::get("derived_and_curated"))
  project_base = GAMBLR.helpers::check_config_value(config::get("project_base"))
  manual_files = dir(paste0(project_base, path_to_files), pattern = ".tsv")

  original_cols <- names(sample_table)
  joined <- sample_table
  for(f in manual_files){
    full = paste0(project_base, path_to_files, f)
    this_data = suppressMessages(read_tsv(full, comment = "#"))
    #TO DO: fix this so it will join on biopsy_id or sample_id depending on which one is present, Done?
    joined = left_join(joined, this_data)
  }

  new_cols <- setdiff(names(joined), setdiff(original_cols, "sample_id"))
  result <- dplyr::select(joined, all_of(new_cols))

  return(result)
}

#' @title Collate Curated SV Results.
#'
#' @description Collate all SV calls from the genome data and summarize for main oncogenes of interest per sample.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return The sample table with additional columns.
#'
#' @import readr dplyr GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' gambl_results_derived = collate_curated_sv_results(gambl_results_derived)
#' }
collate_curated_sv_results = function(sample_table,
                                      seq_type_filter = "genome"){

  new_cols <- compute_curated_sv_results_core(sample_table)
  sample_table = left_join(sample_table, new_cols, by = "sample_id")

  return(sample_table)
}
