#' @title Compute DLBClass labels for a sample scope.
#'
#' @description Core computation behind `collate_dlbclass()`, extracted so
#' it can also be called directly by [collate_results_db()] for just the
#' subset of samples missing from its cache table. Note the join is via
#' `patient_id`, not `sample_id` directly -- samples sharing a patient_id
#' (e.g. multiple biopsies) get the same `dlbclass` value -- but the
#' result is still one row per `sample_id`, matching `sample_table`'s own
#' grain.
#'
#' @param sample_table A data frame with `sample_id` and `patient_id`
#' columns, scoping which samples to compute for.
#'
#' @return A data frame with `sample_id`, `dlbclass`.
#'
#' @import dplyr readr GAMBLR.helpers
#'
#' @keywords internal
#' @noRd
compute_dlbclass_core <- function(sample_table){

    #get paths
    base <- GAMBLR.helpers::check_config_value(
        config::get("project_base")
    )

    dlbclass_s_table <- "icgc_dart/exome_data/dlbclass-s02.csv"

    full_path <- paste0(base, dlbclass_s_table)

    dlbclass <- read_tsv(full_path) %>%
        dplyr::select(
            dlbclass_patient_id = MatchID,
            dlbclass = PredictedCluster
        ) %>%
        dplyr::mutate(
            dlbclass_patient_id = toupper(dlbclass_patient_id),
            dlbclass_patient_id = gsub("_NULLPAIR", "", dlbclass_patient_id),
            dlbclass_patient_id = gsub("-", "_", dlbclass_patient_id)
        )

    result <- dplyr::select(sample_table, sample_id, patient_id) %>%
        dplyr::mutate(
            dlbclass_patient_id = toupper(patient_id),
            dlbclass_patient_id = gsub("-", "_", dlbclass_patient_id)
        ) %>%
        left_join(dlbclass, by = "dlbclass_patient_id") %>%
        dplyr::select(sample_id, dlbclass)

    return(result)
}

#' @title Collate labels from DLBCLass classifier.
#'
#' @description Expand a metadata table horizontally with labels from
#'     Supplemental Table 02 of the DLBCLass paper (PMID: 39680847). This will
#'     only populate these columns for samples where DLBCLass outputs are
#'     availble (Schmitz and Chapuy cohorts), otherwise NA will be returned.
#'
#' @details This is an internal function called by
#'  [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table df with sample ids in the first column. The output of
#'  [GAMBLR.results::get_gambl_metadata] is expected.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr readr GAMBLR.helpers
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#'   sample_table <- get_gambl_metadata(seq_type_filter = "capture")
#'   sample_table <- collate_dlbclass(sample_table = sample_table)
#' }
collate_dlbclass = function(
    sample_table
    ){

    new_cols <- compute_dlbclass_core(sample_table)
    sample_table <- left_join(sample_table, new_cols, by = "sample_id")

    return(sample_table)
}
