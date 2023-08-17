#' @title Collate Derived Results.
#'
#' @description Extract derived results stored in the database (these are usually slower to derive on the fly).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations. Default is FALSE, TRUE is not yet implemented.
#'
#' @return Data frame with one row per sample. Contains the contents of the derived_data table in the database.
#'
#' @import dplyr DBI RMariaDB GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' gambl_results_derived = collate_derived_results(samples_df)
#'
collate_derived_results = function(sample_table,
                                   seq_type_filter = "genome",
                                   from_flatfile = FALSE){

  if(from_flatfile){
    message("not implemented YET")
  }else{
    database_name = GAMBLR.helpers::check_config_value(config::get("database_name"))
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    derived_tbl = dplyr::tbl(con, "derived_data") %>%
      as.data.frame()
  }
  derived_tbl = derived_tbl %>%
    dplyr::select(where( ~!all(is.na(.x)))) %>%
    as.data.frame() #drop the columns that are completely empty

  print(derived_tbl)
  sample_table = dplyr::left_join(sample_table, derived_tbl)
  print(sample_table)
  return(sample_table)
}
