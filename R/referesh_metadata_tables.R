#' @title Refresh Metadata Tables
#'
#' @description Refresh the contents of a metadata table.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package usage.
#'
#' @return Table.
#'
#' @import RMariaDB DBI dplyr
#'
#' @noRd
#'
#' @examples
#' ref_meta = referesh_metadata_tables()
#'
#' @keywords internal
referesh_metadata_tables = function(){

  con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  all_metadata_info = sanity_check_metadata()
  tables = pull(all_metadata_info, table)
  files = pull(all_metadata_info, file)
  for(i in c(1:length(files))){
    refresh_full_table(tables[i], con, files[i])
  }
}
