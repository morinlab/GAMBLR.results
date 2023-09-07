#' @title Refresh Full Table
#'
#' @description Refresh the contents of a database table.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::referesh_metadata_tables], not meant for out-of-package usage.
#'
#' @param table_name Name of table to refresh.
#' @param connection Database connection object.
#' @param file_path Path to the table contents to populate.
#'
#' @return A table.
#'
#' @import DBI RMariaDB readr
#'
#' @noRd
#'
#' @examples
#' refresh_full_table(table_x, con,file_x)
#'
refresh_full_table = function(table_name,
                              connection,
                              file_path){

  table_data = suppressMessages(read_tsv(file_path))
  dbWriteTable(con, table_name, table_data, overwrite = TRUE)
  print(paste("POPULATING table:", table_name, "USING path:", file_path))
}
