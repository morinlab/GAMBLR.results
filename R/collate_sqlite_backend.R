#' @title Connect to the shared collated-results SQLite database.
#'
#' @description Returns a cached DBI connection to the SQLite database that
#' backs [GAMBLR.results::collate_results_db], creating the file (and its
#' parent directory) on first connection if it doesn't exist yet. Unlike
#' `gambl_mutations_db()` in `GAMBLR.data`, this database is writable and
#' local to this repo -- it isn't distributed via GitHub Releases, since
#' its contents (collated results) are internal/pre-publication.
#'
#' @param db_path Optional explicit path to the .db file. Defaults to
#' `results_merged$collated_db` under `project_base`, resolved via
#' `config.yml`.
#'
#' @return A DBIConnection.
#'
#' @import DBI RSQLite GAMBLR.helpers
#' @keywords internal
#' @noRd
gambl_collated_db <- function(db_path = NULL) {
  if (is.null(db_path)) {
    base <- GAMBLR.helpers::check_config_value(config::get("project_base"))
    rel_path <- GAMBLR.helpers::check_config_value(config::get("results_merged")$collated_db)
    db_path <- paste0(base, rel_path)
  }
  cached <- getOption("gamblr.collated.con")
  if (!is.null(cached) && DBI::dbIsValid(cached) &&
      identical(normalizePath(cached@dbname), normalizePath(db_path, mustWork = FALSE))) {
    return(cached)
  }
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  options(gamblr.collated.con = con)
  con
}

#' @title Map an R column's class to a SQLite column type affinity.
#'
#' @description Used only when a column needs to be added to an existing
#' table via `ALTER TABLE ... ADD COLUMN`, since that requires an explicit
#' type (unlike `DBI::dbWriteTable()`'s own type inference on initial table
#' creation).
#'
#' @param x A vector (one column's worth of data).
#'
#' @return A single SQLite type affinity string.
#'
#' @keywords internal
#' @noRd
sqlite_type_affinity <- function(x) {
  if (is.integer(x) || is.logical(x)) return("INTEGER")
  if (is.numeric(x)) return("REAL")
  "TEXT"
}

#' @title Get the set of keys already present in a collate results table.
#'
#' @description Returns the distinct `key_cols` combinations already stored
#' in `table_name`, or a zero-row placeholder if the table doesn't exist
#' yet (i.e. nothing is "already present" -- everything is missing).
#'
#' @param con A DBIConnection, from `gambl_collated_db()`.
#' @param table_name Name of the table to check.
#' @param key_cols Columns identifying a unique result. Default
#' `c("sample_id", "seq_type")`.
#'
#' @return A data frame with just `key_cols`.
#'
#' @import DBI dplyr
#' @keywords internal
#' @noRd
get_existing_collate_keys <- function(con, table_name, key_cols = c("sample_id", "seq_type")) {
  if (!DBI::dbExistsTable(con, table_name)) {
    empty <- as.data.frame(matrix(character(0), ncol = length(key_cols)))
    names(empty) <- key_cols
    return(empty)
  }
  cols_sql <- paste(DBI::dbQuoteIdentifier(con, key_cols), collapse = ", ")
  query <- paste0("SELECT DISTINCT ", cols_sql, " FROM ", DBI::dbQuoteIdentifier(con, table_name))
  DBI::dbGetQuery(con, query)
}

#' @title Write (create or upsert into) a collate results table.
#'
#' @description Creates `table_name` if it doesn't exist yet (schema
#' inferred from `new_data`). If it does exist, adds any columns present in
#' `new_data` but missing from the table (`ALTER TABLE ... ADD COLUMN`,
#' since a `collate_*` helper's output can gain columns over time), then
#' replaces any existing rows sharing a `key_cols` combination with
#' `new_data` and inserts the rest -- delete-matching-rows-then-insert,
#' rather than a raw `INSERT ... ON CONFLICT`, so this doesn't depend on a
#' unique index existing on `key_cols`. Wrapped in a transaction so a
#' failure partway through leaves the table in its prior state rather than
#' a half-written one.
#'
#' @param con A DBIConnection, from `gambl_collated_db()`.
#' @param table_name Name of the table to write to.
#' @param new_data A data frame including `key_cols` plus whatever result
#' columns this `collate_*` helper produced. Must contain at least one row.
#' @param key_cols Columns identifying a unique result. Default
#' `c("sample_id", "seq_type")`.
#'
#' @return Invisibly, `new_data`.
#'
#' @import DBI dplyr
#' @keywords internal
#' @noRd
write_collate_table <- function(con, table_name, new_data, key_cols = c("sample_id", "seq_type")) {
  stopifnot(nrow(new_data) > 0, all(key_cols %in% names(new_data)))

  if (!DBI::dbExistsTable(con, table_name)) {
    DBI::dbWriteTable(con, table_name, new_data)
    return(invisible(new_data))
  }

  existing_cols <- DBI::dbListFields(con, table_name)
  new_cols <- setdiff(names(new_data), existing_cols)
  for (col in new_cols) {
    ddl <- paste0(
      "ALTER TABLE ", DBI::dbQuoteIdentifier(con, table_name),
      " ADD COLUMN ", DBI::dbQuoteIdentifier(con, col),
      " ", sqlite_type_affinity(new_data[[col]])
    )
    DBI::dbExecute(con, ddl)
  }

  DBI::dbBegin(con)
  tryCatch({
    key_combos <- unique(new_data[key_cols])
    where_clauses <- apply(key_combos, 1, function(row) {
      conditions <- mapply(function(col, val) {
        paste0(DBI::dbQuoteIdentifier(con, col), " = ", DBI::dbQuoteString(con, as.character(val)))
      }, key_cols, row)
      paste0("(", paste(conditions, collapse = " AND "), ")")
    })
    delete_sql <- paste0(
      "DELETE FROM ", DBI::dbQuoteIdentifier(con, table_name),
      " WHERE ", paste(where_clauses, collapse = " OR ")
    )
    DBI::dbExecute(con, delete_sql)
    DBI::dbWriteTable(con, table_name, new_data, append = TRUE)
    DBI::dbCommit(con)
  }, error = function(e) {
    DBI::dbRollback(con)
    stop(e)
  })

  invisible(new_data)
}
