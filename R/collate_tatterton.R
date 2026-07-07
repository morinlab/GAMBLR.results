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
}