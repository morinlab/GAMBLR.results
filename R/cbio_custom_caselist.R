#' @title Custom cBioPortal case list.
#'
#' @description Create and a custom case list for easy data subset in cBioPortal.
#'
#' @details Convenience function for specifying custom case lists that can be browsed on cBioPortal.
#' This function takes a set of sample IDs `these_sample_ids` and intersect the IDs with what's available in the study-specific clinical file.
#' This function also extracts the project name for the specified study, i.e the project name that is defined withing the folder specified under the `dir` parameter.
#'
#' @param these_sample_ids A vector of sample IDs to be subset into a case list. Required parameter.
#' @param caselist_name Name of the generated case list (name does not include the file format, this will be added automagically). This parameter is required.
#' @param caselist_description A verbose description of the created case list. Required.
#' @param return_missing_samples Boolean parameter. Set to TRUE to return all sample IDs that are in the desired case list, but not represented in the study specific clinical file. Default is FALSE.
#' @param dir The directory where all study specific files live.
#'
#' @return Nothing.
#' 
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' #get some sample IDs
#' my_samples = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL") %>%
#'  dplyr::filter(cohort == "FL_GenomeCanada") %>%
#'  pull(sample_id)
#'
#' #create case list with selected sample IDs
#' cbio_custom_caselist(these_sample_ids = my_samples,
#'                      caselist_name = "FL_Canada",
#'                      caselist_description = "Follicular Lymphoma from the Genome Canada Study",
#'                      dir = "../path/to/study_directory/")
#' }
#' 
#' @keywords internal
cbio_custom_caselist = function(these_sample_ids,
                                caselist_name,
                                caselist_description,
                                return_missing_samples = FALSE,
                                dir){
  
  #get path to the clinical file holding all sample IDs
  clinical_file = data.table::fread(file = paste0(dir, "data_clinical_samples.txt"), sep = "\t", header = FALSE, skip = 5)
  
  #get project name
  meta_study = data.table::fread(file = paste0(dir, "meta_study.txt"), fill = TRUE)
  project_name = pull(meta_study[2,2])
  
  #pull sample IDs
  clinical_ids = clinical_file %>%
    pull(V2)
  
  #intersect sample IDs from the clinical file with specified sample IDs (these_sample_ids) to ensure no IDs are described in the new case list, that are not represented in the clinical file.
  ids = intersect(clinical_ids, these_sample_ids)
  
  #print how many sample, if any, were removed in the process.
  not_in_clin = setdiff(these_sample_ids, clinical_ids)
  if(!return_missing_samples){
    print(paste0(length(not_in_clin), " Samples not described in the clinical file (data.clinical_samples.txt) for the selected study. To see what samples are not represented, set return_missing_samples = TRUE."))
  }else{
    print(paste0(length(not_in_clin), " Samples not described in the clinical file (data.clinical_samples.txt) for the selected study. The missing samples are: "))
    print(not_in_clin)
  }
  
  #create new case list file
  new_caselist = paste0(dir, "case_lists/", caselist_name, ".txt")
  
  #get sample ID subset belonging to the defined case list
  tabseplist = paste(unique(ids), collapse = "\t")
  
  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, caselist_name),
                   paste0("case_list_name: ", caselist_name),
                   paste0("case_list_description: ", caselist_description),
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))
  
  cat(caselistdata, sep = "\n", file = new_caselist)
}
