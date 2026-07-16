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
#' @export
#'
#' @references Tatterton DJ, Newby ML, Allen JD, et al. The origin, diagnosis,
#'  and prognosis of oligomannose-type diffuse large B-cell lymphoma. Blood.
#'  2025;146(23):2808-2820.
#'
#' @examples
#' \dontrun{
#'   sample_table <- get_gambl_metadata()
#'   sample_table <- collate_tatterton(sample_table = sample_table)
#' }
collate_tatterton = function(
    sample_table
    ){

    # Get Tatterton supplemental table 1b
    base <- GAMBLR.helpers::check_config_value(
        config::get("project_base")
    )

    s1a_path <- paste0(base, "icgc_dart/exome_data/tatterton_s1a.tsv")
    s1b_path <- paste0(base, "icgc_dart/exome_data/tatterton_s1b.tsv")
    assay_path <- paste0(base, "icgc_dart/exome_data/assay_sequencing_data_for_tatterton_s1a.tsv")

    # Rename non-redundant Tatterton columns
    keep_rename <- c(
        Donor_Name = "Donor_Name",
        igseqr_external_AGS_count = "AGS",
        igseqr_external_AGS_location = "AGS_location_type",
        igseqr_external_AGS_motif = "AGS_motif",
        igseqr_external_AGS_codon = "AGS_codon_location",
        igseqr_external_IGHV = "IGHV",
        igseqr_external_IGHJ = "IGHJ",
        igseqr_external_IGHD = "IGHD",
        igseqr_external_IGHC = "IGHC",
        igseqr_external_IGHV_homology_pct = "IGHV_homol_pct",
        lymphgen_tatterton = "LymphGen_call",
        bcl2_tr_tatterton = "BCL2_TR"
    )

    # Read Tatterton supplemental tables and keep only columns in keep_rename
    s1a <- suppressMessages(readr::read_tsv(s1a_path, col_types = cols(.default = "c"))) %>%
        dplyr::select(dplyr::any_of(unname(keep_rename)))
    s1b <- suppressMessages(readr::read_tsv(s1b_path, col_types = cols(.default = "c"))) %>%
        dplyr::select(dplyr::any_of(unname(keep_rename)))

    tatterton <- dplyr::bind_rows(s1a, s1b) %>%
        dplyr::distinct(Donor_Name, .keep_all = TRUE)

    # Read in Assay_Sequencing table excerpt to support S1a sample annotation
    assay_dlc <- suppressMessages(read_tsv(assay_path, col_types = cols(.default = "c")))

    # Create Assay_Sequencing map for linking S1a samples
    map_s1a <- assay_dlc %>%
        dplyr::filter(!grepl("mi_RNA|miRNA", Sample_ID, ignore.case = TRUE)) %>% # Filter out miRNA data
        dplyr::distinct(patient_id, biopsy_id, Donor_Name) %>%
        dplyr::rename(tatterton_rna_biopsy = biopsy_id) %>%
        dplyr::filter(Donor_Name %in% s1a$Donor_Name)
    
    map_s1b <- s1b %>%
        dplyr::transmute(patient_id = Donor_Name, Donor_Name,
                         tatterton_rna_biopsy = NA_character_)
    
    id_map <- dplyr::bind_rows(map_s1a, map_s1b) %>%
        dplyr::distinct(patient_id, tatterton_rna_biopsy, .keep_all = TRUE)

    tatterton_calls <- id_map %>%
        dplyr::inner_join(tatterton, by = "Donor_Name") %>%
        dplyr::distinct(patient_id, tatterton_rna_biopsy, .keep_all = TRUE)
    
    tatterton_patients <- unique(tatterton_calls$patient_id)

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