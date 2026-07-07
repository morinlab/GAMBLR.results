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
}