#' @title Compute QC metrics for a sample scope.
#'
#' @description Core computation behind `collate_qc_results()`, extracted
#' so it can also be called directly by [collate_results_db()] for just
#' the subset of samples missing from its cache table. Reads one QC file
#' per (`unix_group`, `seq_type`) combination present in `sample_table` --
#' bounded by how many distinct combinations are in scope, not by sample
#' count, so unlike SSM this doesn't have a "loads everything regardless
#' of scope" problem.
#'
#' @param sample_table A data frame with `sample_id`, `seq_type`, and
#' `unix_group` columns, scoping which samples to compute for.
#'
#' @return A data frame with `sample_id` plus whatever QC metric columns
#' the source files contain.
#'
#' @import dplyr readr glue purrr GAMBLR.helpers
#'
#' @keywords internal
#' @noRd
compute_qc_results_core <- function(sample_table){

    # Fool-proof against unexpected seq types
    supported_seq_types <- sample_table %>%
        filter(
            seq_type %in% c("genome", "capture")
        )

    # get paths
    base <- GAMBLR.helpers::check_config_value(config::get("project_base"))
    qc_template <- GAMBLR.helpers::check_config_value(config::get("qc_met"))

    paths <- expand.grid(
        unix_group = unique(supported_seq_types$unix_group),
        seq_type_filter = unique(supported_seq_types$seq_type),
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
        read_tsv(
            qc_path_full,
            show_col_types = FALSE
        )%>%
        mutate(
            unix_group = unix_group
        ) %>%
        dplyr::rename(
            sample_id = UID,
            seq_type = SeqType
        )

    }
    )

    # unix_group/seq_type are metadata already tracked on the caller's own
    # sample_table -- drop them here so this function only contributes its
    # own genuinely new QC columns, not a duplicate copy of metadata.
    qc_metrics <- dplyr::select(qc_metrics, -unix_group, -seq_type)

    # every requested sample_id gets a row back, even with NA QC values if
    # it wasn't found in the relevant file(s).
    result <- dplyr::select(sample_table, sample_id) %>%
        left_join(qc_metrics, by = "sample_id")

    return(result)
}

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

    new_cols <- compute_qc_results_core(sample_table)
    sample_table = left_join(sample_table, new_cols, by = "sample_id")

    #print n samples with QC metrics (at least one non-NA metric column)
    qc_cols <- setdiff(names(new_cols), "sample_id")
    n_with_qc <- sum(apply(dplyr::select(new_cols, all_of(qc_cols)), 1, function(row) any(!is.na(row))))
    message(paste("QC metrics for", n_with_qc, "samples retrieved."))

    return(sample_table)
}
