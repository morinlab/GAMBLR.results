#' @title Collate oligomannose type DLBCL annotations from Tatterton et al 2025.
#'
#' @description Expand a sample table horizontally with the Mann-type
#'      (oligomannose-type) DLBCL classifications from Tatterton et al., Blood 2025
#'      (supplemental Table 1B). These columns are populated for samples that have
#'      Tatterton data (the dlbcl_schmitz cohort) and NA otherwise.
#'      A case is Mann-type when it has BOTH an acquired N-glycosylation
#'      site (AGS) in the CDR AND a follicular-lymphoma (FL) signature (EZB subtype or a BCL2
#'      translocation). AGS calls are patient-level (RNA-seq derived), so the
#'      Tatterton table (keyed on Donor_Name) is left-joined onto the sample table
#'      `patient_id`. Samples absent from the Tatterton table (outside dlbcl_schmitz cohort)
#'      receive NA for the added columns. `manntype` and `FL_signature` are POS/NEG factors.
#'
#' @details This is an internal function called by
#'  [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table df with sample ids in the first column. The output of
#'  [GAMBLR.results::get_gambl_metadata] is expected. The column patient_id is also required to join the information from Tatterton table S1b.
#'
#' @return The sample table with `sample_id`, `manntype`, `FL_signature`, `AGS`,
#'  `AGS_location_type`, and `AGS_Motif` appended, followed by the remaining metadata 
#'  and Tatterton table columns.
#'
#' @import dplyr readr GAMBLR.helpers
#'
#' @noRd
#'
#' @references Tatterton DJ, Newby ML, Allen JD, et al. The origin, diagnosis,
#'  and prognosis of oligomannose-type diffuse large B-cell lymphoma. Blood.
#'  2025;146(23):2808-2820.
#'
#' @examples
#' \dontrun{
#'   sample_table <- get_gambl_metadata(seq_type_filter = "capture")
#'   sample_table <- collate_tatterton(sample_table = sample_table)
#' }
collate_tatterton = function(
    sample_table
    ){

    # Get Tatterton supplemental table 1b
    base <- GAMBLR.helpers::check_config_value(
        config::get("project_base")
    )

    tatterton_table <- "icgc_dart/exome_data/tatterton_s1b.csv"

    full_path <- paste0(base, tatterton_table)

    # Read Tatterton s1b and derive the Mann-type calls per patient
    # FL_signature: POS if LymphGen_call contains "EZB" or BCL2_TR == "POS"
    # manntype: POS if AGS in CDR (AGS >= 1 and AGS_location_type contains "CDR") AND FL_signature == "POS"
    # All other Tatterton columns are kept and pushed to the far right using everything()
    tatterton <- suppressMessages(read_csv(full_path)) %>%
        dplyr::rename(patient_id = Donor_Name) %>% # Rename Donor_Name in tatterton to patient_id so it's shared with sample_table
        dplyr::mutate(
            FL_signature = if_else( # Create FL_signature column. POS if LymphGen_call from tatterton contains "EZB" OR BCL2_TR == "POS"
                grepl("EZB", toupper(coalesce(LymphGen_call, ""))) |
                    toupper(coalesce(BCL2_TR, "")) == "POS",
                "POS", "NEG"
            ),
            manntype = if_else( # Create manntype column. POS if AGS >=1 AND AGS is in the CDR AND FL_signature == "POS"
                grepl("CDR", toupper(coalesce(AGS_location_type, ""))) &
                    coalesce(suppressWarnings(as.numeric(AGS)), 0) >= 1 &
                    FL_signature == "POS",
                "POS", "NEG"
            ),
            across(c(FL_signature, manntype), ~ factor(.x, levels = c("NEG", "POS")))
        ) %>%
        dplyr::distinct(patient_id, .keep_all = TRUE)

    # Preserve the incoming column names so the original and Tatterton columns can be reordered after the join
    sample_table_cols <- colnames(sample_table)
 
    # Left join onto the sample table (Donor_Name links to patient_id). Samples not in the Tatterton table get NA for the added columns
    sample_table <- left_join(sample_table, tatterton, by = "patient_id")
 
    # Reorder columns to sample_id, Mann-type annotations, the remaining original sample_table columns, and the remaining Tatterton columns
    first_cols <- c(
        "sample_id", "manntype", "FL_signature",
        "AGS", "AGS_location_type", "AGS_Motif"
    )
    sample_table <- sample_table %>%
        dplyr::select(
            dplyr::any_of(first_cols),
            dplyr::any_of(setdiff(sample_table_cols, first_cols)),
            dplyr::everything()
        )

    return(sample_table)
}