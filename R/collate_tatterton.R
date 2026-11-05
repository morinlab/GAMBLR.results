#' @title Collate oligomannose-type (Mann-type) DLBCL annotations from Tatterton et al. 2025.
#'
#' @description Expand a sample table horizontally with Mann-type (oligomannose-type)
#'      DLBCL classifications from Tatterton et al., Blood 2025, covering BOTH
#'      supplemental cohorts: Table S1A (BCCA) and Table S1B (NCI). A case is
#'      Mann-type when it has BOTH an acquired N-glycosylation site (AGS) in the CDR
#'      AND a follicular-lymphoma (FL) signature (EZB subtype or a BCL2 translocation).
#'
#'      The AGS call (CDR N-glycosylation site) is IgSeqR-derived from mRNA and has no
#'      GAMBL equivalent, so it is taken from Tatterton and retained under
#'      `igseqr_external_*` column names. Two Mann-type derivations are returned:
#'      `manntype` / `FL_signature` derive the FL signature from GAMBL's own `lymphgen`
#'      and `bcl2_ba` (GAMBL as source of truth), while `manntype_tatterton` /
#'      `FL_signature_tatterton` reproduce the original published derivation using
#'      Tatterton's `LymphGen_call` and `BCL2_TR` (retained as `lymphgen_tatterton`
#'      and `bcl2_tr_tatterton`). Both share the same Tatterton AGS half.
#'
#'      Tatterton's AGS is anchored to the specific mRNA biopsy it was derived from.
#'      S1A (BCCA) donors are `DLC####` and are linked to that biopsy via an
#'      Assay_Sequencing crosswalk (microRNA libraries excluded); a sample is
#'      annotated only if it is on that mRNA biopsy (its mRNA plus any DNA sharing it),
#'      so a patient's other biopsies are NOT annotated. S1B (NCI) donors equal
#'      `patient_id` and are linked patient-level, since their RNA and DNA share
#'      biopsies. Samples outside both cohorts receive NA for the added columns.
#'      `manntype`, `FL_signature`, `manntype_tatterton`, and `FL_signature_tatterton`
#'      are POS/NEG factors.
#'
#' @details This is an internal function called by
#'  [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table df with `sample_id` in the first column, expected to be the
#'  output of [GAMBLR.results::get_gambl_metadata] with NO seq_type filter applied
#'  (mRNA samples are required to anchor the biopsy link). The columns `patient_id`,
#'  `biopsy_id`, and `seq_type` are required.
#' @param compare_derivations logical; if TRUE, print how many `manntype` and
#'  `FL_signature` calls differ between the GAMBL-derived and Tatterton-derived
#'  versions, and list the flipped samples. Default FALSE.
#'
#' @return The sample table with Mann-type annotations appended
#'  (`manntype`, `FL_signature`, `manntype_tatterton`, `FL_signature_tatterton`),
#'  the IgSeqR-derived AGS columns (`igseqr_external_AGS_count`,
#'  `igseqr_external_AGS_location`, `igseqr_external_AGS_motif`,
#'  `igseqr_external_AGS_codon`), a `tatterton_link` provenance column, followed by
#'  the remaining metadata and Tatterton columns. Samples outside the Tatterton
#'  cohorts receive NA for the added columns.
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

    # Create df with which samples from sample table will be annotated with Tatterton columns
    # S1a is linked using tatterton_rna_biopsy = biopsy_id
    # S1b is linked using Donor_Name = patient_id
    matched <- sample_table %>%
        dplyr::filter(patient_id %in% tatterton_patients) %>%
        dplyr::left_join(
            dplyr::select(tatterton_calls, patient_id, tatterton_rna_biopsy),
            by = "patient_id", relationship = "many-to-many"
        ) %>%
        dplyr::group_by(sample_id, seq_type, biopsy_id, patient_id) %>%
        dplyr::summarise(
            on_tatt_biopsy = any(biopsy_id == tatterton_rna_biopsy, na.rm = TRUE),
            no_anchor = all(is.na(tatterton_rna_biopsy)),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            tatterton_link = dplyr::case_when(
                on_tatt_biopsy ~ "tatterton_biopsy",
                no_anchor ~ "s1b_patient_link",
                TRUE ~ NA_character_
            )
        ) %>%
        dplyr::filter(!is.na(tatterton_link)) %>%
        dplyr::select(sample_id, seq_type, biopsy_id, tatterton_link)
    
    sample_table_cols <- colnames(sample_table)

    # Attach Tatterton columns to matching rows from sample table
    tatterton_join <- tatterton_calls %>%
        dplyr::select(-tatterton_rna_biopsy) %>%
        dplyr::rename(dplyr::any_of(keep_rename)) %>%
        dplyr::distinct(patient_id, .keep_all = TRUE)

    sample_table <- sample_table %>%
        dplyr::left_join(matched, by = c("sample_id", "seq_type", "biopsy_id")) %>%
        dplyr::left_join(tatterton_join, by = "patient_id", na_matches = "never") %>%
        dplyr::mutate(
            dplyr::across(
                dplyr::any_of(names(keep_rename)[-1]),   # incl. lymphgen/bcl2 tatterton
                ~ dplyr::if_else(is.na(tatterton_link), NA, .x)
            )
        )
    
    # Annotate manntype and FL_signature
    # FL_signature = POS if LymphGen subtype is EZB or BCL2 translocation is present
    # manntype = POS if FL_signature == "POS" and AGS count >= 1 and AGS location == "CDR"
    sample_table <- sample_table %>%
        dplyr::mutate(
            # GAMBL-derived manntype annotations
            FL_signature = dplyr::if_else(
                !is.na(tatterton_link) & (
                    grepl("EZB", toupper(dplyr::coalesce(lymphgen, ""))) |
                    toupper(dplyr::coalesce(bcl2_ba, "")) == "POS"),
                "POS", "NEG"
            ),
            manntype = dplyr::if_else(
                !is.na(tatterton_link) &
                    grepl("CDR", toupper(dplyr::coalesce(igseqr_external_AGS_location, ""))) &
                    dplyr::coalesce(suppressWarnings(as.numeric(igseqr_external_AGS_count)), 0) >= 1 &
                    FL_signature == "POS",
                "POS", "NEG"
            ),
            # Tatterton-derived manntype annotations (original)
            FL_signature_tatterton = dplyr::if_else(
                !is.na(tatterton_link) & (
                    grepl("EZB", toupper(dplyr::coalesce(lymphgen_tatterton, ""))) |
                    toupper(dplyr::coalesce(bcl2_tr_tatterton, "")) == "POS"),
                "POS", "NEG"
            ),
            manntype_tatterton = dplyr::if_else(
                !is.na(tatterton_link) &
                    grepl("CDR", toupper(dplyr::coalesce(igseqr_external_AGS_location, ""))) &
                    dplyr::coalesce(suppressWarnings(as.numeric(igseqr_external_AGS_count)), 0) >= 1 &
                    FL_signature_tatterton == "POS",
                "POS", "NEG"
            ),
            # Set manntype annotations to NA for rows without a match in Tatterton data
            dplyr::across(
                c(FL_signature, manntype, FL_signature_tatterton, manntype_tatterton),
                ~ dplyr::if_else(is.na(tatterton_link), NA_character_, .x)
            ),
            # Set manntype annotations to character
            dplyr::across(
                c(FL_signature, manntype, FL_signature_tatterton, manntype_tatterton),
                ~ factor(.x, levels = c("NEG", "POS"))
            )
        )
 
    # Reorder columns
    first_cols <- c(
        "sample_id", "Donor_Name", "manntype", "FL_signature", "manntype_tatterton",
        "FL_signature_tatterton", "igseqr_external_AGS_count", "igseqr_external_AGS_location",
        "igseqr_external_AGS_motif", "igseqr_external_AGS_codon"
    )

    sample_table <- sample_table %>%
        dplyr::select(
            dplyr::any_of(first_cols),
            dplyr::any_of(setdiff(sample_table_cols, first_cols)),
            dplyr::everything()
        )

    return(sample_table)
}
