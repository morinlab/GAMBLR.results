#' @title Collate Cellularity/Purity from Battenberg results.
#'
#' @description Expand a metadata table horizontally with cellularity and ploidy
#'  estimated from battenberg. This will only populate these columns for samples
#'  where battenberg outputs are availble, otherwise NA will be returned.
#'
#' @details This is an internal function called by
#'  [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table df with sample ids in the first column. The output of
#'  [GAMBLR.results::get_gambl_metadata] is expected.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr readr glue GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' sample_table <- get_gambl_metadata()
#' bat <- collate_battenberg_purity(sample_table = sample_table)
#'
collate_battenberg_purity = function(
    sample_table
    ){

    #get paths
    base <- GAMBLR.helpers::check_config_value(
        config::get("project_base")
    )
    bat_template <- GAMBLR.helpers::check_config_value(
        config::get("results_flatfiles")$cnv$battenberg_cellularity
    )

    bat_path <- sample_table %>%
        rowwise() %>%
        mutate(
            path_to_bat = paste0(
                base,
                glue::glue(bat_template)
            )
        )

    ploidy <- apply(
        bat_path %>% select(path_to_bat, sample_id), 1, function(row) {
            file_path <- row[1]
            if (file.exists(file_path)) {
                suppressMessages(read_tsv(file_path)) %>%
                    mutate(sample_id = row[2]) %>%
                    as.data.frame()
            } else {
                return(NULL)
            }
    })

    ploidy <- do.call(bind_rows, ploidy)

    sample_table <- left_join(sample_table, ploidy)

    return(sample_table)
}
