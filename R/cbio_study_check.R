#' @title Study Check (cBioPortal).
#'
#' @description Helper function for checking integrity of study files.
#'
#' @details This function was designed to ensure that all the sample IDs described in the maf are actually present in the clinical files.
#' If this is not the case, the function will notify the user what samples are found in the case list that are not described in the clinical file.
#' The function then sub-sets the case list to only include samples from the clinical file.
#' Note that the `project_name` has to match what is specified for the previously run functions (i.e [GAMBLR.results::cbio_setup_study], [GAMBLR.results::cbio_setup_fusions] and [GAMBLR.results::cbio_finalize_study]).
#'
#' @param data_clinical_samples_path Path to clinical file.
#' @param data_fusions_path Path to data_fusion file from setup_fusions.
#' @param cases_fusions_path Path to cases_fusion from setup_fusions.
#' @param cases_all_path Path to cases_all from setup_study.
#' @param cases_sequenced_path Path to cases_sequenced from setup_study.
#' @param project_name Project name, should match what is specified under setup_study/setup_fusions.
#' @param out_dir Directory with all study related files, the only argument that needs to be specified, given that paths to all generated study files are not changed from default.
#'
#' @return Nothing.
#'
#' @import dplyr readr
#' @export
#'
#' @examples
#' \dontrun{
#' samples_not_in_clinical = cbio_study_check(out_dir = "GAMBLR/cBioPortal/instance01/")
#' }
cbio_study_check = function(data_clinical_samples_path = "data_clinical_samples.txt",
                            data_fusions_path = "data_fusions.txt",
                            cases_fusions_path = "case_lists/cases_fusion.txt",
                            cases_all_path = "case_lists/cases_all.txt",
                            cases_sequenced_path = "case_lists/cases_sequenced.txt",
                            project_name = "gambl_genome",
                            out_dir){

  #read clinical file (skip header)
  data_clinical_samples = data.table::fread(file = paste0(out_dir, data_clinical_samples_path), sep = "\t", header = FALSE, skip = 5)

  #cases fusions
  cases_fusion = data.table::fread(file = paste0(out_dir, cases_fusions_path), sep = "	", skip = 4, header = FALSE)

  #transform data and strip irrelevant characters
  cases_fusion = t(cases_fusion) %>%
    as.data.frame() %>%
    mutate_at("V1", str_replace, "case_list_ids: ", "")

  #return samples that are present in case_lists but not in clinical file
  not_in_case_lists = setdiff(cases_fusion$V1, data_clinical_samples$V2)

  if(length(not_in_case_lists) > 0){
    #print message
    message("Warning! Samples found in case lists that are not described in the clinical file (data_clincaal_samples.txt). The following samples will be removed from all case lists and data fusions:")
    print(not_in_case_lists)

    #intersect with clinical file to remove sample ids that are not represented (in clinical file)
    cases_fusions_check = intersect(cases_fusion$V1, data_clinical_samples$V2)

    #put updated sample ids back into case list and overwrite original file
    caselist_fusions = paste0(out_dir, cases_fusions_path)

    tabseplist_fusions = paste(cases_fusions_check, collapse = "\t")

    caselistdata_fusions = c(paste0("cancer_study_identifier: ", project_name),
                             paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                             paste0(c("case_list_ids:", tabseplist_fusions),collapse = " "))

    cat(caselistdata_fusions, sep = "\n", file = paste0(out_dir, cases_fusions_path))

    #data fusions
    data_fusions = data.table::fread(file = paste0(out_dir, data_fusions_path), sep = "\t", header = TRUE)

    #remove all samples not in clinical file that are present in data fusions
    data_fusions = data_fusions[!grepl(paste(not_in_case_lists, collapse = "|"), data_fusions$Tumor_Sample_Barcode),]

    #overwrite data_fusions.txt
    data_fusions_out = paste0(out_dir, data_fusions_path)
    write_tsv(data_fusions, data_fusions_out)

    #cases all
    cases_all = data.table::fread(file = paste0(out_dir, cases_sequenced_path), sep = "	", skip = 4, header = FALSE)

    #transform data and strip irrelevant characters
    cases_all = t(cases_all) %>%
      as.data.frame() %>%
      mutate_at("V1", str_replace, "case_list_ids: ", "")

    #intersect with clinical file to remove sample ids that are not represented (in clinical file)
    cases_all_check = intersect(cases_all$V1, data_clinical_samples$V2)

    #put updated sample ids back into case list and overwrite original file
    caselist_all = paste0(out_dir, cases_all_path)

    tabseplist_all = paste(cases_all_check, collapse = "\t")

    caselistdata_all = c(paste0("cancer_study_identifier: ", project_name),
                         paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                         paste0(c("case_list_ids:", tabseplist_all),collapse = " "))

    cat(caselistdata_all, sep = "\n", file = paste0(out_dir, cases_all_path))

    #cases sequenced
    cases_sequenced = data.table::fread(file = paste0(out_dir, cases_sequenced_path), sep = "	", skip = 4, header = FALSE)

    #transform data and strip irrelevant characters
    cases_sequenced = t(cases_sequenced) %>%
      as.data.frame() %>%
      mutate_at("V1", str_replace, "case_list_ids: ", "")

    #intersect with clinical file to remove sample ids that are not represented (in clinical file)
    cases_sequenced_check = intersect(cases_sequenced$V1, data_clinical_samples$V2)

    #put updated sample ids back into case list and overwrite original file
    caselist_sequenced = paste0(out_dir, cases_sequenced_path)

    tabseplist_sequenced = paste(cases_sequenced_check, collapse = "\t")

    caselistdata_sequenced = c(paste0("cancer_study_identifier: ", project_name),
                               paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                               paste0(c("case_list_ids:", tabseplist_sequenced),collapse = " "))

    cat(caselistdata_sequenced, sep = "\n", file = paste0(out_dir, cases_sequenced_path))
  }else{
    message("No additional samples were found in case lists that are not described in the clinical file (data_clinical_samples.txt)")
  }
}
