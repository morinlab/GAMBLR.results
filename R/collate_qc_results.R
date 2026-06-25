#' @title Collate Quality Control Results.
#'
#' @description Expand a metadata table horizontally with quality control metrics.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table df with sample ids in the first column.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr readr glue purrr GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#'   qc_metrics = collate_qc_results(sample_table = sample_table)
#' }
collate_qc_results = function(sample_table){

    #get paths
    base <- GAMBLR.helpers::check_config_value(config::get("project_base"))
    qc_template <- GAMBLR.helpers::check_config_value(config::get("qc_met"))

    paths <- expand.grid(
        unix_group = unique(sample_table$unix_group),
        seq_type_filter = unique(sample_table$seq_type),
        stringsAsFactors = FALSE
    ) %>%
    mutate(
        qc_path = pmap_chr(
            list(unix_group, seq_type_filter),
            ~ glue(
                qc_template,
                unix_group = ..1,
                seq_type_filter = ..2
            )
        ),
        qc_path_full = paste0(base, qc_path)
    )

    qc_metrics <- pmap_dfr(
        paths,
    \(unix_group, seq_type_filter, qc_path, qc_path_full) {
        read_tsv(qc_path_full) %>%
        mutate(
            unix_group = unix_group,
            seq_type = seq_type_filter
        )
    }
    )
    # read in qc data, rename sample id and seq type columns
    qc_metrics <- pmap_dfr(
        paths,
    \(unix_group, seq_type_filter, qc_path, qc_path_full) {
        read_tsv(qc_path_full) %>%
        mutate(
            unix_group = unix_group
        ) %>%
        dplyr::rename(
            sample_id = UID,
            seq_type = SeqType
        )

    }
    )


    # join sample table and QC data
    sample_table = left_join(sample_table, qc_metrics)

    #print n samples with QC metrics
    qc_samples = length(unique(qc_metrics$sample_id))
    message(paste("QC metrics for", qc_samples, "samples retrieved."))

    return(sample_table)
}
