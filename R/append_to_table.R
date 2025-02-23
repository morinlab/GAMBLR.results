#' @title Append To Table.
#'
#' @description Housekeeping function to add results to a table.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package usage.
#'
#' @param table_name The name of the database table to update/populate.
#' @param data_df A dataframe of values to load into the table.
#'
#' @return A table.
#'
#' @import RMariaDB DBI GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#'   table_up = append_to_table("my_table", "my_df")
#' }
append_to_table = function(table_name,
                           data_df){

  db = GAMBLR.helpers::check_config_value(config::get("database_name"))
  con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  dbWriteTable(con, table_name, table_data, append = TRUE)
}
