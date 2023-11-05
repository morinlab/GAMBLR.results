#' @title Finalize Study (cBioPortal).
#'
#' @description Finish setting up a new cBioPortal instance or updating an existing portal data set.
#'
#' @details This function should be run as the last (or third step) in setting up a new cBioPortal instance.
#' The functions that should be run prior to these functions are; [GAMBLR.results::setup_study] and [GAMBLR.results::setup_fusions].
#' [GAMBLR.results::finalize_study] creates all the necessary tables and metadata files and case lists that are required to import a new study into cBioPortal.
#' Note, that all parameter arguments used in this function have to match the same parameter arguments for the previously run functions (`setup_study` and `setup_fusions`).
#' This function allows the user to specify additional fields from the collated metadata file (besides the "standard" fields).
#' For more information on how to use, see `metacols` and related parameters (`metacol_names`, `metacol_types`, and `meta_prior`).
#'
#' @param seq_type_filter the seq type you are setting up a study for, default is "genome".
#' @param short_name A concise name for your portal project.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param description A verbose description of your data set.
#' @param cancer_type Cancer types included in study, default is "mixed".
#' @param these_sample_ids A vector of all the sample_id that were included in any of the data files for cBioPortal (i.e the output from `setup_study` and `setup_fusions`).
#' @param overwrite Flag to specify that files should be overwritten if they exist. Default is TRUE.
#' @param metacols Optional, specify any additional metadata/collate_result fields to be included in the cBioPortal metadata. Extra columns are specified as a vector of characters. If not provided, only "standard" metadata fields will be kept.
#' @param metacol_names Separately specify the names of the columns given to `metacols` as a vector of characters. The number of elements needs to match the total number of columns specified with `metacols`. Required parameter if `metacols` is being called.
#' @param metacol_types Specify the data type for selected metadata columns as a vector of characters. The number of elements specified needs to match the number of selected columns with `metacols`. Acceptable values are; STRING, NUMBER and BOOLEAN. Required parameter if `metacols` is being called.
#' @param metacol_prior Explicitly state the priority of selected metadata columns as a vector of characters. A higher number indicates a higher priority. Required parameter if `metacols` is being called.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return Nothing.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' #basic usage
#' cbio_finalize_study(these_sample_ids = c(ids, fusion_ids), out_dir = "GAMBLR/cBioPortal/instance01/")
#'
#' #advanced usage
#' #get some samples
#' all_meta = get_gambl_metadata()
#' meta_sub = head(all_meta, 5)
#' my_samples = pull(meta_sub, sample_id)
#'
#' #create a clinical file with additional collated metadata fields
#' cbio_finalize_study(these_sample_ids = my_samples,
#'                     out_dir = "../",
#'                     metacols = c("MeanCorrectedCoverage", "total_ssm"),
#'                     metacol_names = c("Mean Corrected Coverage", "Total SSM"),
#'                     metacol_types = c("NUMBER", "NUMBER"),
#'                     metacol_prior = c("2", "1"))
#' }
#'
cbio_finalize_study = function(seq_type_filter = "genome",
                               short_name = "GAMBL",
                               human_friendly_name = "GAMBL data",
                               project_name = "gambl_genome",
                               description = "GAMBL data from genome",
                               cancer_type = "mixed",
                               these_sample_ids,
                               overwrite = TRUE,
                               metacols,
                               metacol_names,
                               metacol_types,
                               metacol_prior,
                               out_dir){

  #define standard columns
  these_columns = c("patient_id", "sample_id", "pathology",
                    "EBV_status_inf", "cohort", "time_point",
                    "ffpe_or_frozen", "myc_ba", "bcl6_ba",
                    "bcl2_ba", "COO_consensus", "DHITsig_consensus", "lymphgen")

  these_names = c("Patient Identifier", "Sample Identifier", "Subtype",
                  "EBV status", "Cohort", "Time point",
                  "FFPE", "MYC_BA", "BCL6_BA",
                  "BCL2_BA", "COO", "DHITsig", "LymphGen")

  these_types = c(rep("STRING", 13))

  these_priorities = c("1", "1", "3",
                       "2", "4", "2",
                       "2", "2", "2",
                       "2", "2", "2", "4")

  #add any extra columns, if such are specified.
  if(!missing(metacols)){
    metacols = c(these_columns, metacols)
  }else{
    metacols = these_columns
  }

  if(!missing(metacol_names)){
    metacol_names = c(these_names, metacol_names)
  }else{
    metacol_names = these_names
  }

  if(!missing(metacol_types)){
    metacol_types = c(these_types, metacol_types)
  }else{
    metacol_types = these_types
  }

  if(!missing(metacol_prior)){
    metacol_prior = c(these_priorities, metacol_prior)
  }else{
    metacol_prior = these_priorities
  }

  #create necessary files
  #create case list
  caselist = paste0(out_dir, "case_lists/cases_sequenced.txt")

  tabseplist = paste(unique(these_sample_ids), collapse = "\t")

  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist),collapse = " "))

  cat(caselistdata, sep = "\n", file = caselist)

  #create case list all
  caselist_all = paste0(out_dir, "case_lists/cases_all.txt")

  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ",project_name, "_allcases"), "case_list_name: Samples sequenced", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))

  cat(caselistdata, sep = "\n", file = caselist_all)

  #meta samples
  #prepare and write out the relevant metadata
  clinsamp = paste0(out_dir, "data_clinical_samples.txt")
  meta_samples = collate_results(seq_type_filter = seq_type_filter, join_with_full_metadata = TRUE) %>%
    dplyr::filter(sample_id %in% these_sample_ids) %>%
    dplyr::select(all_of(metacols))

  colnames(meta_samples) = toupper(colnames(meta_samples))

  #wrangle to get header based on what columns are selected
  #deal with column names
  if(length(metacols) == length(metacol_names)){
    col_names_tmp = paste(metacol_names, collapse = '\t')
    col_names = paste0("#", col_names_tmp)
  }else{
    message("The number of column names specified with `metacol_names` does not match the number of columns specified with `metacols`")
    message(paste0(length(metacol_names), " elements are specified with `metacol_names`, and there are ", length(metacols), " elements in the `metacols` parameter"))
    stop()
  }

  #deal with column data type
  if(length(metacols) == length(metacol_types)){
    col_types_tmp = paste(metacol_types, collapse = '\t')
    col_types = paste0("#", col_types_tmp)
  }else{
    message("The number of elements in `metacol_types` does not match the number of columns specified with `metacols`:")
    message(paste0(length(metacol_types), " elements are specified with `metacol_types`, and there are ", length(metacols), " elements in the `metacols` parameter"))
    stop()
  }

  #deal with column priority
  if(length(metacols) == length(metacol_prior)){
    col_prior_tmp = paste(metacol_prior, collapse = '\t')
    col_prior = paste0("#", col_prior_tmp)
  }else{
    message("The number of elements in `metacol_prior` does not match the number of columns specified with `metacols`")
    message(paste0(length(metacol_prior), " elements are specified with `metacol_prior`, and there are ", length(metacols), " elements in the `metacols` parameter"))
    stop()
  }

  #construct the header
  header = paste0(col_names, "\n", col_names, "\n", col_types, "\n", col_prior, "\n")

  cat(header, file = clinsamp)

  write.table(meta_samples, file = clinsamp, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)

  #create clinical meta data
  #first get the patient_id list
  clinpat = paste0(out_dir, "data_clinical_patient.txt")
  patient_ids = pull(meta_samples, PATIENT_ID)

  all_outcomes = get_gambl_outcomes(time_unit = "month", censor_cbioportal = TRUE, patient_ids = patient_ids, complete_missing = TRUE) %>%
    dplyr::select(c("patient_id", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS", "age", "sex"))

  colnames(all_outcomes) = toupper(colnames(all_outcomes))

  header = paste0("#Patient Identifier\tOverall Survival Status\tOverall Survival (Months)\tDisease Free Status\tDisease Free (Months)\tAGE\tSEX\n",
                  "#Patient Identifier\tOverall Survival Status\tOverall Survival (Months)\tDisease Free Status\tDisease Free (Months)\tAge\tSex\n",
                  "#STRING\tSTRING\tNUMBER\tSTRING\tNUMBER\tNUMBER\tSTRING\n",
                  "#1\t1\t1\t1\t1\t1\t1\n")

  cat(header, file = clinpat)

  write.table(all_outcomes, file = clinpat, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
}
