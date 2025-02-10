#' @title Find Expected Outputs.
#'
#' @description Get the details including file paths for the anticipated outputs from a pipeline or tool.
#'
#' @details This function takes a tool or pipeline with `tool_name` and the unix group with `unix_group` and returns information such as paths to individual files.
#' Optionally, the user can provide an already loaded data frame with all the file details (`targ_df`).
#' for more information and examples, refer to the parameter descriptions as well as function examples.
#'
#' @param targ_df Optionally provide a data frame with all file details.
#' @param tool_name The tool or pipeline that generated the files (should be the same for all). Acceptable values are manta and gridss.
#' @param unix_group The unix group (should be the same for all).
#' @param filename_end_pattern Optionally specify a pattern to search for the files among a longer set of files in the outputs.
#' @param update_db Set to TRUE to overwrite any existing rows in the table for this tool/unix_group combination.
#' @param target_path Path to targets.
#'
#' @import dplyr readr RMariaDB stringr DBI tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' #get paths to unmatched manta bedpe files
#'ex_outs = find_expected_outputs(tool_name = "manta",
#'                                unix_group = "gambl",
#'                                filename_end_pattern = "unmatched.somaticSV.bedpe")
#'
#' @keywords internal
find_expected_outputs = function(targ_df,
                                 tool_name,
                                 unix_group,
                                 filename_end_pattern,
                                 update_db = FALSE,
                                 target_path){

  repo_base = GAMBLR.helpers::check_config_value(config::get("repo_base"))
  if(missing(target_path)){
    target_path = paste0(repo_base, "targets/", tool_name, "--", unix_group)
  }

  if(tool_name == "manta"){
    if(missing(targ_df)){
      filename_end_pattern = ".somaticSV.bedpe"

      targ_df = suppressMessages(read_tsv(target_path, col_names = c("file"))) %>%
        dplyr::filter(str_detect(file, pattern = filename_end_pattern))

      targ_df = mutate(targ_df, file_path = paste0(repo_base, file)) %>%
        separate(file, sep = "/", into = c("results", "unix_group", "tool_version", "outputs", "type", "seq_genome", "detail", "filename"))

      targ_df = separate(targ_df, seq_genome, sep = "--", into = c("seq_type", "genome_build")) %>%
        dplyr::select(-results, -type, -detail, -outputs) %>%
        separate(filename, sep = "--", into = c("tumour_sample_id", "normal_sample_id", "pairing_status")) %>%
        mutate(pairing_status = str_remove(pairing_status, filename_end_pattern)) %>%
        separate(tool_version, sep = "-", into = c("tool_name", "tool_version"))
    }
    targ_df = targ_df %>%
      mutate(file_timestamp = file.info(file_path)$mtime)

    targ_df$output_type = "bedpe"
    print(targ_df)
    print(target_path)

  }else if(tool_name == "gridss"){
    filename_end_pattern = ".gridss_somatic_filtered.bedpe"

    targ_df = suppressMessages(read_tsv(target_path, col_names = c("file"))) %>%
      dplyr::filter(str_detect(file, pattern = filename_end_pattern))

    targ_df = mutate(targ_df, file_path = paste0(repo_base, file)) %>%
      separate(file, sep = "/", into = c("results", "unix_group", "tool_version", "outputs", "type", "seq_genome", "detail", "filename"))

    targ_df = separate(targ_df, seq_genome, sep = "--", into = c("seq_type", "genome_build")) %>%
      dplyr::select(-results, -type, -detail, -outputs) %>%
      separate(filename, sep = "--", into = c("tumour_sample_id", "normal_sample_id", "pairing_status")) %>%
      mutate(pairing_status = str_remove(pairing_status, filename_end_pattern)) %>%
      separate(tool_version, sep = "-", into = c("tool_name", "tool_version")) %>%
      mutate(file_timestamp = file.info(file_path)$mtime)

    targ_df$output_type = "bedpe"
  }
  if(update_db){
    database_name = GAMBLR.helpers::check_config_value(config::get("database_name"))
    con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    table_name = GAMBLR.helpers::check_config_value(config::get("tables")$files)
    message(paste("updating", table_name,"in", database_name))
    #clear all files for this tool/unix_group combination
    update_q = paste0("DELETE from ", table_name, " WHERE tool_name = \"", tool_name, "\" and unix_group = \"", unix_group, "\" ;")
    print(update_q)
    something = dbReadTable(conn = con, table_name)

    summarized = something %>%
      group_by(unix_group) %>%
      tally()

    #merge with incoming data
    print(summarized)
    #to_write = rbind(something, targ_df)
    dbExecute(con, update_q, immediate = TRUE)
    dbWriteTable(con, table_name, targ_df, append = TRUE)

    #dbExecute(con, update_q, immediate = TRUE)
    #dbWriteTable(con, table_name, expected_manta, append = TRUE)
    DBI::dbDisconnect(con)
  }
  return(targ_df)
}
