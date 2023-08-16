#' @title Get GAMBL Outcomes.
#'
#' @description Get the patient-centric clinical metadata.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::get_gambl_metadata], not meant for out-of-package usage.
#'
#' @param patient_ids Vector of patient IDs.
#' @param time_unit Return follow-up times in one of three time units: year, month or day. Default is "year".
#' @param censor_cbioportal Optionally request the censoring to be encoded in the specific style required by cBioPortal. Default is FALSE.
#' @param complete_missing Optionally fill in any gaps to ensure we have values for every patient (censor at 0 if missing). Default is FALSE.
#' @param from_flatfile Optionally set to FALSE to use the database to get the survival data. Default is TRUE.
#'
#' @return Data frame with one row for each patient_id.
#'
#' @import tidyr dplyr readr RMariaDB DBI
#'
#' @noRd
#'
#' @examples
#' outcome_df = get_gambl_outcomes()
#'
#' @export
get_gambl_outcomes = function(patient_ids,
                              time_unit = "year",
                              censor_cbioportal = FALSE,
                              complete_missing = FALSE,
                              from_flatfile = TRUE){

  if(from_flatfile){
    outcome_flatfile = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("table_flatfiles")$outcomes))

    #check for missingness
    if(!file.exists(outcome_flatfile)){
      print(paste("missing: ", outcome_flatfile))
      message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
    }

    all_outcome = suppressMessages(read_tsv(outcome_flatfile))

  }else{
    db = check_config_value(config::get("database_name"))
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
    all_outcome = dplyr::tbl(con, "outcome_metadata") %>%
      as.data.frame()

    DBI::dbDisconnect(con)
  }
  if(!missing(patient_ids)){
    all_outcome = all_outcome %>%
      dplyr::filter(patient_id %in% patient_ids)

    if(complete_missing){
      #add NA values and censored outcomes for all missing patient_ids
      all_outcome = all_outcome %>%
        complete(patient_id = patient_ids, fill = list(OS_YEARS = 0, PFS_years = 0, TTP_YEARS = 0, DSS_YEARS = 0, CODE_OS = 0, CODE_PFS = 0, CODE_DSS = 0, CODE_TTP = 0))
    }
  }
  if(time_unit == "month"){
    all_outcome = all_outcome %>%
      mutate(OS_MONTHS = OS_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(PFS_MONTHS = PFS_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(TTP_MONTHS = TTP_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(DSS_MONTHS = DSS_YEARS * 12)

    all_outcome = all_outcome %>%
      dplyr::select(-c("OS_YEARS", "PFS_YEARS", "TTP_YEARS", "DSS_YEARS"))

  }else if(time_unit == "day"){
    all_outcome = all_outcome %>%
      mutate(OS_DAYS = OS_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(PFS_DAYS = PFS_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(TTP_DAYS = TTP_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(DSS_DAYS = DSS_YEARS * 365)

    all_outcome = all_outcome %>%
      dplyr::select(-c("OS_YEARS", "PFS_YEARS", "TTP_YEARS", "DSS_YEARS"))
  }

  #if necessary, convert the censoring into the cBioPortal format for OS and PFS
  if(censor_cbioportal){
    all_outcome$OS_STATUS = as.character(all_outcome$CODE_OS)
    all_outcome = all_outcome %>%
      mutate(OS_STATUS = case_when(OS_STATUS == "0" ~ "0:LIVING", OS_STATUS == "1"~"1:DECEASED"))

    all_outcome$DFS_STATUS = as.character(all_outcome$CODE_PFS)
    all_outcome = all_outcome %>%
      mutate(DFS_STATUS = case_when(DFS_STATUS == "0" ~ "0:DiseaseFree", DFS_STATUS == "1"~"1:Recurred/Progressed"))

    all_outcome = all_outcome %>%
      mutate(all_outcome, DFS_MONTHS = PFS_MONTHS)
  }
  all_outcome = all_outcome %>%
    mutate(is_adult = ifelse(age < 20, "Pediatric", "Adult"))

  return(all_outcome)
}
