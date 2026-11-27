#' @title Get GAMBL Outcomes.
#'
#' @description Get the patient-centric clinical metadata.
#'
#' @details One row per patient. To add these columns to sample metadata,
#' join by `patient_id` only (`by = "patient_id"`): without it, dplyr joins on
#' every shared column (e.g. `sex`), which silently fails to match rows whose
#' values differ between the two tables.
#'
#' Cleaning applied to the outcome table:
#' \itemize{
#'   \item `sex` is "M", "F" or NA (sources spell it e.g. "Male"/"Female").
#'   \item `IPI` and `FLIPI` codes of -1 (not assessable) are NA.
#'   \item `tx_primary_class` collapses the primary treatment (`tx_primary`)
#'     into a few classes: the first systemic regimen decides (surgery,
#'     splenectomy or radiation given before it are skipped) and rituximab
#'     counts when a course of it is given (maintenance or consolidation
#'     alone don't). Classes: "R-CHOP-like" (CHOP, CEOP, ECHOP, COEP, ACOP or
#'     BPVACOP with rituximab), "BR" (bendamustine + rituximab), "R-CVP",
#'     "R-mono" (rituximab alone), "observation", "radiation" (radiation
#'     only), "chemo, no rituximab" (the CHOP-like or CVP regimens without
#'     rituximab) and "other" (anything else, including intensive regimens
#'     such as CODOX-M/IVAC or DA-EPOCH, fludarabine-based or targeted
#'     therapy, and free-text entries like "Systemic therapy"; free-text
#'     "Radiation therapy" alone is "radiation"; "R-CHOP-like", as recorded
#'     for published cohorts, is "R-CHOP-like" and "non-R-CHOP-like" is
#'     "other"). NA when unknown.
#'   \item `tx_ritu_maintenance`: TRUE if the primary treatment included
#'     rituximab maintenance, FALSE if a coded regimen didn't, NA when the
#'     treatment isn't a coded regimen.
#' }
#'
#' @param patient_ids Vector of patient IDs.
#' @param time_unit Return follow-up times in one of three time units: year, month or day. Default is "year".
#' @param censor_cbioportal Optionally request the censoring to be encoded in the specific style required by cBioPortal. Default is FALSE.
#' @param complete_missing Optionally add a row for every requested patient missing from the outcome table, censored at time 0. Patients already in the table keep their missing values as NA. Default is FALSE.
#' @param from_flatfile Optionally set to FALSE to use the database to get the survival data. Default is TRUE.
#'
#' @return Data frame with one row for each patient_id.
#'
#' @import tidyr dplyr readr DBI GAMBLR.helpers
#'
#' @examples
#' \dontrun{
#'   outcome_df = get_gambl_outcomes()
#'   # join by patient_id only
#'   metadata_with_outcomes = dplyr::left_join(
#'                             get_gambl_metadata(),
#'                             dplyr::select(outcome_df, -sex),
#'                             by = "patient_id")
#' }
#' @export
get_gambl_outcomes = function(patient_ids,
                              time_unit = "year",
                              censor_cbioportal = FALSE,
                              complete_missing = FALSE,
                              from_flatfile = TRUE){

  if(from_flatfile){
    outcome_flatfile = paste0(check_config_and_value("repo_base"),
                              check_config_and_value("table_flatfiles$outcomes"))

    #check for missingness
    if(!file.exists(outcome_flatfile)){
      print(paste("missing: ", outcome_flatfile))
      message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
    }

    all_outcome = suppressMessages(read_tsv(outcome_flatfile))

  }
  all_outcome = clean_gambl_outcomes(all_outcome)

  if(!missing(patient_ids)){
    all_outcome = all_outcome %>%
      dplyr::filter(patient_id %in% patient_ids)

    if(complete_missing){
      #add rows with censored outcomes (time 0) for patient_ids missing from the table; explicit = FALSE keeps
      #the NAs of patients already in it (e.g. no PFS recorded) as NA
      all_outcome = all_outcome %>%
        complete(patient_id = patient_ids, fill = list(OS_YEARS = 0, PFS_YEARS = 0, TTP_YEARS = 0, DSS_YEARS = 0,
                                                       TTT_YEARS = 0, CODE_OS = 0, CODE_PFS = 0, CODE_DSS = 0,
                                                       CODE_TTP = 0, CODE_TTT = 0),
                 explicit = FALSE)
    }
  }
  all_outcome = convert_outcome_time_unit(all_outcome, time_unit)

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

# Follow-up times (the *_YEARS columns) in months or days instead of years (365.25 days per year)
convert_outcome_time_unit = function(outcome, time_unit = "year"){
  if(!time_unit %in% c("year", "month", "day")){
    stop("time_unit must be one of \"year\", \"month\" or \"day\"")
  }
  if(time_unit == "year") return(outcome)
  year_cols = intersect(c("OS_YEARS", "PFS_YEARS", "TTP_YEARS", "DSS_YEARS", "TTT_YEARS"), colnames(outcome))
  factor = if(time_unit == "month") 12 else 365.25
  suffix = if(time_unit == "month") "_MONTHS" else "_DAYS"
  for(col in year_cols){
    outcome[[sub("_YEARS$", suffix, col)]] = outcome[[col]] * factor
  }
  dplyr::select(outcome, -dplyr::all_of(year_cols))
}

# Cleaning shared with GAMBLR.data's bundled outcomes: sex as M/F, IPI/FLIPI -1 as NA, coarse treatment classes
clean_gambl_outcomes = function(outcome){
  if("sex" %in% colnames(outcome)) outcome$sex = standardize_sex(outcome$sex)
  for(col in intersect(c("IPI", "FLIPI"), colnames(outcome))){
    v = suppressWarnings(as.integer(outcome[[col]]))
    v[v %in% -1] = NA
    outcome[[col]] = v
  }
  if("tx_primary" %in% colnames(outcome)){
    outcome$tx_primary_class = classify_primary_treatment(outcome$tx_primary)
    outcome$tx_ritu_maintenance = primary_treatment_ritu_maintenance(outcome$tx_primary)
  }
  outcome
}

# "M", "F" or NA from the spellings used across sources (M/F, Male/Female, male/female)
standardize_sex = function(sex){
  s = tolower(trimws(as.character(sex)))
  dplyr::case_when(s %in% c("m", "male") ~ "M",
                   s %in% c("f", "female") ~ "F",
                   TRUE ~ NA_character_)
}

# The tokens of a coded regimen (e.g. "SURG-CHOP3-RITU-RAD"): upper case, a trailing note like " (PHOENIX)"
# dropped, without course counts ("CHOP3" -> "CHOP")
treatment_tokens = function(tx){
  tx = sub(" \\(.*\\)$", "", toupper(trimws(as.character(tx))))
  lapply(strsplit(tx, "-", fixed = TRUE), function(t) sub("[0-9]+$", "", trimws(t)))
}

# TRUE for regimens written in the BC Cancer code style (e.g. "CHOP6-RITU"), FALSE for free text
# ("Systemic therapy, Radiation therapy"), "UNKN" and NA
is_coded_treatment = function(tx){
  s = trimws(as.character(tx))
  !is.na(s) & grepl("^[A-Z0-9/-]+( \\(.*\\))?$", s) & !toupper(s) %in% "UNKN"
}

classify_primary_treatment = function(tx){
  chop_like = c("CHOP", "CEOP", "ECHOP", "COEP", "ACOP", "BPVACOP")
  cvp_like = c("CVP", "CVPP", "CVPR")
  radiation = c("RAD", "ISRT", "IRRT")
  skip = c("SURG", "SPLE", radiation)
  out = rep("other", length(tx))
  out[is.na(tx) | toupper(trimws(as.character(tx))) %in% "UNKN"] = NA_character_
  # entries that aren't coded regimens: free text from MOHCCN (only radiation alone can be classed) and the
  # R-CHOP-like yes/no of published cohorts (e.g. Chapuy et al.)
  free = tolower(trimws(as.character(tx)))
  out[free %in% "radiation therapy"] = "radiation"
  out[free %in% "r-chop-like"] = "R-CHOP-like"
  coded = is_coded_treatment(tx)
  tokens = treatment_tokens(tx)
  for(i in which(coded)){
    t = tokens[[i]]
    # a course of rituximab: a RITU token, or RITU in a combined token like "IBR/PLB" + "RITU"
    rituximab = any(t == "RITU") || any(grepl("(^|/)RITU($|/)", t))
    systemic = t[!t %in% c(skip, "OBSE")]
    first = if(length(systemic)) systemic[1] else NA_character_
    out[i] = if(is.na(first)){
      if(any(t == "OBSE")) "observation" else if(any(t %in% radiation)) "radiation" else "other"
    }else if(first %in% chop_like){
      if(rituximab) "R-CHOP-like" else "chemo, no rituximab"
    }else if(first == "BEND"){
      if(rituximab) "BR" else "other"
    }else if(first %in% cvp_like){
      if(rituximab || first == "CVPR") "R-CVP" else "chemo, no rituximab"
    }else if(first == "RITU" && all(systemic %in% c("RITU", "RITUMAIN"))){
      "R-mono"
    }else{
      "other"
    }
  }
  out
}

primary_treatment_ritu_maintenance = function(tx){
  out = rep(NA, length(tx))
  coded = is_coded_treatment(tx)
  tokens = treatment_tokens(tx)
  out[coded] = vapply(tokens[coded], function(t) any(t == "RITUMAIN"), logical(1))
  out
}
