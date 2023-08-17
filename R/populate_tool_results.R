#' @title Populate Tool Results.
#'
#' @description Populate the database with the per-sample summarized results of various tools.
#'
#' @details this function is still in draft mode, export to NAMESPACE has been removed for now.
#'
#' @param tool_name Name of the tool to get the results for.
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' results = populate_tool_results(tool_name = "slims_3")
#' }
#' 
populate_tool_results = function(tool_name){

  #IMPORTANT TODO: This function should only ever work with samples that exist in the metadata
  # Perhaps it should report any excluded outputs in case they need to be deleted from the main output directories
  matched_analyses = unlist(GAMBLR.helpers::check_config_value(config::get("analyses")$matched))
  print(matched_analyses)
  database_name = GAMBLR.helpers::check_config_value(config::get("database_name"))
  genome_builds = unlist(strsplit(GAMBLR.helpers::check_config_value(config::get("genome_builds")), ","))
  groups = unlist(strsplit(GAMBLR.helpers::check_config_value(config::get("unix_groups")), ","))
  for(analysis_type in names(matched_analyses)){
    tool_name = matched_analyses[analysis_type]
    message(paste("populating results for", tool_name))
    populate_each_tool_result(tool = tool_name, genome_builds, groups)
  }
}
