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
            dlbclass_patient_id = gsub("-", "_", dlbclass_patient_id),
            dlbclass = paste0("C", dlbclass)
        )


    sample_table <- left_join(
        sample_table %>%
            mutate(
                dlbclass_patient_id = toupper(patient_id),
                dlbclass_patient_id = gsub("-", "_", dlbclass_patient_id)
            ),
        dlbclass
    ) %>%
        dplyr::select(-dlbclass_patient_id)

    return(sample_table)
}
