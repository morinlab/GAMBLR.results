#' @title Fetch Output Files.
#'
#' @description Get output files from a set of conditions.
#'
#' @details This function lets the user specify multiple conditions for returning result subsets.
#' First, specify the name of the tool with `tool`, then set the seq type (`this_seq_type`) to either genome or capture,
#' together with the genome build (`genome_build`). A data frame will be returned with one row per file and sample IDs together with GAMBL wildcards.
#'
#' @param tool Name of tool.
#' @param unix_group The unix group of the sample set.
#' @param base_path Either the full or relative path to where all the results directories are for the tool e.g. "gambl/sequenza_current".
#' @param results_dir Directory with results.
#' @param this_seq_type Either genome or capture.
#' @param build Default is hg38.
#' @param search_pattern File-extensions search pattern.
#'
#' @return A data frame with one row per file and sample IDs parsed from the file name along with other GAMBL wildcards.
#'
#' @import dplyr tibble tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' ex_outs = fetch_output_files(tool = "manta",
#'                              base_path = "gambl/sequenza_current",
#'                              this_seq_type = "capture",
#'                              build = "hg38")
#' }
#' @keywords internal
fetch_output_files = function(tool,
                              unix_group,
                              base_path,
                              results_dir = "99-outputs",
                              this_seq_type = "genome",
                              build = "hg38",
                              search_pattern = "cellularity_ploidy.txt"){

  if(!grepl("^/", base_path)){
    project_base = GAMBLR.helpers::check_config_value(config::get("project_base",config="default"))
    local_project_base = GAMBLR.helpers::check_config_value(config::get("project_base"))
    #project_base = "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
    local_base_path = paste0(local_project_base, base_path)
    base_path = paste0(project_base, base_path)

  }
  if(tool == "battenberg"){
    results_path = paste0(base_path, "/", results_dir, "/seg/", this_seq_type, "--projection/")
    local_results_path = paste0(local_base_path, "/", results_dir, "/seg/", this_seq_type, "--projection/")
  }else{
    results_path = paste0(base_path, "/", results_dir, "/", this_seq_type, "--projection/")
    local_results_path = paste0(local_base_path, "/", results_dir, "/", this_seq_type, "--projection/")
  }

  #This still fails when a matching file isn't found. No clue why this doesn't work
  if(tool == "sequenza"){

    print(paste0("using results in: ", results_path))
    print("THIS CAN BE SLOW!")

    unnested_df = unnested_df  %>%
      mutate(full_path = paste0(results_path, short_path, "/filtered/sequenza_alternative_solutions.txt"))

    named = pull(unnested_df, full_path)

    found_files = tibble(filename = lapply(named, file.exists)) %>%
      unnest_longer(filename)

    new_df = cbind(unnested_df, found_files) %>%
      dplyr::filter(found_files == TRUE)

    print(head(new_df))
  }else if(tool == "battenberg"){
    results_path = paste0(base_path, "/", results_dir, "/seg/", this_seq_type, "--", build,"/")
    all_files = dir(results_path, pattern = search_pattern)
    print(results_path)
    print(search_pattern)
    #extract tumour and normal ID
    all_tumours = unlist(lapply(all_files, function(x){tumour = unlist(strsplit(x, "--"))[1]}))
    all_normals = unlist(lapply(all_files, function(x){tumour = unlist(strsplit(x, "--"))[2]}))
    all_files = unlist(lapply(all_files, function(x){paste0(results_path, x)}))
    new_df = data.frame(tumour_sample_id = all_tumours, normal_sample_id = all_normals, full_path = all_files)
    new_df = mutate(new_df, normal_sample_id = gsub(normal_sample_id, pattern = "_subclones.igv.seg", replacement = ""))
  }else if(tool == "battenberg_ploidy"){
    #results_path = paste0(base_path, "/", results_dir, "/", this_seq_type, "--", build, "/")
    search_pattern = "_cellularity_ploidy.txt"
    print("THIS CAN BE SLOW!")
    unnested_df = unnested_df %>%
      mutate(full_path = paste0(results_path, short_path))

    named = pull(unnested_df, full_path)
    found_files = tibble(filename = lapply(named, function(x){dir(x, pattern = search_pattern)[1]})) %>%
      unnest_longer(filename) #%>% mutate(path = paste0(full_path, filename))

    found_files = cbind(found_files, unnested_df) %>%
      dplyr::filter(!is.na(filename)) %>%
      mutate(full_path = paste0(full_path, "/", filename))

    return(found_files)
  }else if(tool == "manta"){
    tool_results_path = GAMBLR.helpers::check_config_value(config::get("results_directories")$manta)
    search_pattern=".bedpe"
    new_df = find_files_extract_wildcards(tool_name = "manta", genome_build = c("hg38", "grch37"), search_pattern = ".bed")

    #details to include
    #assemble_file_details = function(file_paths, tool_name, unix_group, sample_ids, output_type = "ploidy", is_production = "yes")
  }
  return(new_df)
}
