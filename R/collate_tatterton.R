#' @title Collate oligomannose type DLBCL annotations from Tatterton et al 2025.
#'
#' @description Return the Tatterton et al. (Blood 2025, supplemental Table 1B) NCI
#'      DLBCL cohort samples with their Mann-type (oligomannose-type) classification,
#'      joined to GAMBL metadata. Takes no arguments; always builds from the
#'      dlbcl_schmitz cohort (the only cohort with Tatterton data).
#'
#' @details A case is Mann-type when it has both an acquired N-glycosylation site
#'      (AGS) in the CDR AND a follicular-lymphoma (FL) signature (EZB subtype
#'      or a BCL2 translocation). AGS calls are patient-level (derived from
#'      RNA-seq), so the Tatterton table (keyed on `Donor_Name`) is joined onto
#'      the metadata `patient_id` via an inner join, returning only samples
#'      present in the Tatterton cohort. `manntype` and `FL_signature` are POS/NEG
#'      factors. Note the join is patient-level, so a patient with both a capture
#'      and an mrna sample yields one row per seq_type carrying the same call.
#'      Use `dplyr::distinct(patient_id, manntype)` for per-patient counts.
#'
#' @return A data frame of dlbcl_schmitz samples that are in the Tatterton
#'      cohort, with information from Tatterton supplemental table 1b columns plus
#'      `manntype` and `FL_signature` annotations, reordered with the Mann-type
#'      annotation information up front. One row per sample (not per patient).
#'
#' @import dplyr readr stringr GAMBLR.helpers
#'
#' @references Tatterton DJ, Newby ML, Allen JD, et al. The origin, diagnosis,
#'      and prognosis of oligomannose-type diffuse large B-cell lymphoma. Blood.
#'      2025;146(23):2808-2820.
#'
#' @examples
#' \dontrun{
#'   manntype_meta <- collate_tatterton()
#' }
collate_tatterton = function(){
    # Get Tatterton supplemental table 1b
    base <- GAMBLR.helpers::check_config_value(config::get("project_base"))
    tatterton_file <- paste0(
        base,
        "icgc_dart/exome_data/tatterton_s1b.csv"
    )
    tatterton_full <- suppressMessages(read_csv(tatterton_file, col_names = TRUE))

    # Fetch dlbcl_schmitz cohort metadata
    # Keep samples in Tatterton using inner_join
    joined <- get_gambl_metadata() %>%
        filter(cohort == "dlbcl_schmitz") %>%
        inner_join(tatterton_full, by = c("patient_id" = "Donor_Name")) # Donor_Name in Tatterton s1b links to patient_id in dlbcl_schmitz metadata
    
    # Mutate FL_signature (POS if BCL2_TR == "POS" or LymphGen_call contains "EZB")
    # and manntype (POS if AGS in CDR and FL_signature == "POS").
    joined <- joined %>%
        mutate(
            FL_signature = if_else(
                str_detect(toupper(coalesce(LymphGen_call, "")), "EZB") |
                    toupper(coalesce(BCL2_TR, "")) == "POS",
                "POS", "NEG"
            ),
            manntype = if_else(
                str_detect(toupper(coalesce(AGS_location_type, "")), "CDR") &
                    coalesce(suppressWarnings(as.numeric(AGS)), 0) >= 1 &
                    FL_signature == "POS",
                "POS", "NEG"
            ),
            across(c(FL_signature, manntype), ~ factor(.x, levels = c("NEG", "POS"))) # Convert column types to factor
        )

    # Warn if any Tatterton patient had no matching sample in the metadata
    n_expected <- dplyr::n_distinct(tatterton_full$Donor_Name)
    n <- dplyr::n_distinct(joined$patient_id)
    if(n < n_expected){
        warning(n_expected - n, " Tatterton patient(s) had no matching ",
                "sample in the metadata and are not in the output.")
    }

    # Reorder columns. First, sample_id, manntype data first, and key metadata
    first_cols <- c(
        "sample_id", "manntype", "FL_signature",
        "AGS", "AGS_location_type", "AGS_Motif",
        "patient_id", "pathology", "seq_type",
        "capture_space", "pairing_status", "ffpe_or_frozen",
        "biopsy_id", "genome_build", "COO_Class",
        "LymphGen_call", "BCL2_TR", "MYC_TR",
        "DZsig", "bcl6_ba", "bcl2_ba",
        "myc_ba"
    )


    # The remaining Tatterton columns
    tatterton_cols <- setdiff(
        colnames(tatterton_full),
        c("Donor_Name", first_cols)
    )

    # More relevant GAMBL metadata columns
    end_cols <- c("lymphgen", "Tumor_Sample_Barcode")

    # The rest of the GAMBL metadata columns are in everything()
    joined <- joined %>%
        select(
            any_of(first_cols),
            any_of(tatterton_cols),
            any_of(end_cols),
            everything()
        )

    return(joined)
}