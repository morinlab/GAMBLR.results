#' @title Find Files Extract Wildcards
#'
#' @description Get wildcards for a set of samples.
#'
#' @details Specify the file extension with `search_pattern` and the seq type, unix group, and genome build and the function will return a tibble with sample wildcards.
#'
#' @param tool_results_path Optional parameter, path to results.
#' @param search_pattern Search pattern.
#' @param genome_build Genome projection to be used.
#' @param seq_type Default is genome.
#' @param unix_group Default value is gambl.
#' @param tool_name Name of the tool.
#'
#' @return A tibble with wildcards.
#'
#' @import dplyr tidyr tibble
#' @export
#'
#' @examples
#' file_details_manta = find_files_extract_wildcards(tool_name = "manta",
#'                                                   genome_build = c("hg38", "grch37"),
#'                                                   search_pattern = ".bed")
#'
find_files_extract_wildcards = function(tool_results_path,
                                        search_pattern,
                                        genome_build,
                                        seq_type = "genome",
                                        unix_group = "gambl",
                                        tool_name){

  project_base = check_config_value(config::get("project_base"))
  if(missing(tool_results_path)){
    tool_results_paths = check_config_value(config::get("results_directories"))
    tool_results_path = tool_results_paths[[tool_name]]
  }
  results_paths = paste0(project_base, unix_group, "/", tool_results_path, "genome--", genome_build, "/somaticSV/")

  found_files = tibble(filename = lapply(results_paths, function(x){dir(x, pattern = search_pattern)})) %>%
    mutate(path = results_paths, genome_build = genome_build) %>%
    unnest_longer(filename) %>%
    dplyr::filter(!is.na(filename)) %>%
    mutate(file_path = paste0(path, filename)) %>%
    mutate(pairing_status = case_when(grepl("--unmatched", filename) ~ "unmatched", TRUE ~"matched")) %>%
    dplyr::mutate(sample_id = strsplit(filename, "--")) %>%
    unnest_wider(sample_id, names_sep = "-") %>%
    dplyr::rename(tumour_sample_id = `sample_id-1`, normal_sample_id = `sample_id-2`) %>%
    dplyr::mutate(tool_name = tool_name, seq_type = seq_type, unix_group = unix_group) %>%
    dplyr::select(tumour_sample_id, unix_group, tool_name, seq_type, genome_build, file_path, pairing_status, normal_sample_id)

  return(found_files)
}
