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
#' \dontrun{
#' sample_table <- get_gambl_metadata()
#' bat <- collate_battenberg_purity(sample_table = sample_table)
#' }
#' @keywords internal
collate_battenberg_purity = function(
    sample_table
    ){

    # paths from the config, via get_battenberg (genomes with a matched normal)
    bat_files <- suppressMessages(
        get_battenberg(sample_table, file_types = "cellularity_ploidy")
    )
    bat_files <- bat_files[bat_files$cellularity_ploidy_exists, , drop = FALSE]

    ploidy <- lapply(seq_len(nrow(bat_files)), function(i) {
        suppressMessages(read_tsv(bat_files$cellularity_ploidy_file[i])) %>%
            mutate(sample_id = bat_files$sample_id[i]) %>%
            as.data.frame()
    })

    ploidy <- do.call(bind_rows, ploidy)

    sample_table <- left_join(sample_table, ploidy)

    return(sample_table)
}
