#' @title Assemble File Details.
#'
#' @description Update the database by appending to the gambl_files table.
#'
#' @details Specify the file paths with `file_paths` followed by the name of the tool (`tool_name`).
#' Next, set the output type (e.g seq, maf, etc.) and unix group (should be the same for all).
#' Lastly, specify the sample IDs with `these_sample_ids`.
#' For more information on how to use the optional parameters, refer to the parameter descriptions.
#'
#' @param file_details_df Optionally supply the data frame directly instead (e.g. from [GAMBLR::find_files_extract_wildcards]).
#' @param file_paths A vector of full file paths, e.g. the output of dir.
#' @param tool_name The tool or pipeline that generated the files (should be the same for all).
#' @param unix_group The unix group (should be the same for all).
#' @param these_sample_ids A vector of sample_id the same length and in the same order as the file paths.
#' @param output_type The file type to distinguish different output file types from the same pipeline (e.g. seg, maf, ploidy).
#' @param is_production Boolean parameter. Default is yes.
#'
#' @return Updates the database by appending to the gambl_files table. Use with caution!
#'
#' @import tibble RMariaDB DBI tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' assemble_file_details(file_paths = c(one.maf, another.maf), 
#'                       tool_name = "manta",
#'                       unix_group = "genome",
#'                       output_type = "maf",
#'                       these_sample_ids = c(one_sample, another_sample))
#' }
#' 
assemble_file_details = function(file_details_df,
                                 file_paths,
                                 tool_name,
                                 unix_group,
                                 these_sample_ids,
                                 output_type = "ploidy",
                                 is_production = "yes"){

  #the gambl_files table contains
  #sample_id, unix_group, tool_name, tool_version, seq_type, genome_build, is_production, output_type, is_lifted_over, pairing_status
  # is_production, output_type, file_path, file_timestamp
  if(missing(file_details_df)){
  file_details_df = tibble(sample_id = these_sample_ids, unix_group = unix_group, tool_name = tool_name, output_type = output_type, is_production = is_production,
                           file_timestamp = lapply(file_paths, function(x){as.character(file.mtime(x))}),
                           file_path = file_paths) %>%
                              unnest_longer(file_timestamp)
  }

  #date_info = file.mtime(file_path)
  con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  dbWriteTable(con, "gambl_files", file_details_df, append = TRUE)
}
