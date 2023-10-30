#' @title Setup Study (cBioPortal).
#'
#' @description Initialize a new cBioPortal instance or update existing portal data set, can also be used to retrieve sample ids included in study.
#'
#' @details This function internally calls [GAMBLR.results::get_coding_ssm] to retrieve coding mutations to be included in the study (if `overwrite = TRUE`).
#' In addition, this function also creates and sets up the proper folder hierarchy and writes the files necessary to import a new cBioPortal study.
#' Before a study is ready to be imported to cBioPortal, the user also needs to run [GAMBLR.results::setup_fusions] and [GAMBLR.results::finalize_study].
#' Optionally the user can also run [GAMBLR.results::study_check] to ensure all samples described by the "clinical" file are included in the study.
#' Also, note that the parameters chosen for this function have to match the same parameters called for any subsequent study function calls.
#'
#' @param seq_type_filter the seq type you are setting up a study for, default is "genome".
#' @param short_name A concise name for your portal project.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param description A verbose description of your data set.
#' @param overwrite Flag to specify that files should be overwritten if they exist. Default is TRUE.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return A vector of sample_id for the patients that have been included.
#'
#' @import dplyr readr
#' @export
#'
#' @examples
#' #Setup study and save included ids as a vector of characters:
#' \dontrun{
#' ids = cbio_setup_study(out_dir = "GAMBLR/cBioPortal/instance01/")
#' }
#'
cbio_setup_study = function(seq_type_filter = "genome",
                            short_name = "GAMBL",
                            human_friendly_name = "GAMBL data",
                            project_name = "gambl_genome",
                            description = "GAMBL data from genome",
                            overwrite = TRUE,
                            out_dir){
  
  cancer_type="mixed"

  #set up the new directory
  if (!file.exists(out_dir)){
    dir.create(out_dir)
    dir.create(paste0(out_dir, "case_lists"))
  }else{
    dir.create(paste0(out_dir, "case_lists"))
  }

  #create necessary files
  #meta study
  meta_study = paste0(out_dir, "meta_study.txt")

  meta_study_content = paste0("type_of_cancer: ", cancer_type, "\n",
                              "cancer_study_identifier: ", project_name, "\n",
                              "name: ", human_friendly_name, "\n",
                              "short_name: ", short_name, "\n",
                              "description: ", description, "\n",
                              "add_global_case_list: true\n")

  cat(meta_study_content, file = meta_study)


  #meta clinical samples
  meta_clinical_samples = paste0(out_dir, "meta_clinical_samples.txt")

  meta_clinical_samples_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                         "genetic_alteration_type: CLINICAL\n",
                                         "datatype: SAMPLE_ATTRIBUTES\n",
                                         "data_filename: data_clinical_samples.txt")

  cat(meta_clinical_samples_content, file = meta_clinical_samples)

  #meta clinical patients
  meta_clinical_patients = paste0(out_dir, "meta_clinical_patient.txt")

  meta_clinical_patients_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                          "genetic_alteration_type: CLINICAL\n",
                                          "datatype: PATIENT_ATTRIBUTES\n",
                                          "data_filename: data_clinical_patient.txt")

  cat(meta_clinical_patients_content, file = meta_clinical_patients)

  #meta mutations
  meta_mutations = paste0(out_dir, "meta_mutations_extended.txt")
  data_mutations = "data_mutations_extended.maf"

  meta_mutations_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                  "genetic_alteration_type: MUTATION_EXTENDED\n",
                                  "datatype: MAF\n",
                                  "stable_id: mutations\n",
                                  "show_profile_in_analysis_tab: true\n",
                                  "profile_description: Mutation data\n",
                                  "profile_name: Mutations\n",
                                  "data_filename: ", data_mutations, "\n",
                                  "swissprot_identifier: name\n")

  cat(meta_mutations_content, file = meta_mutations)

  if(overwrite){
    #create the actual MAF file by querying the database using the API
    coding_ssms = get_coding_ssm(this_seq_type = seq_type_filter)
    data_mutations_full = paste0(out_dir, "data_mutations_extended.maf")
    write_tsv(coding_ssms, data_mutations_full, na = "")
  }else{
    #read in the MAF instead
    coding_ssms = data.table::fread(file = data_mutations_full,
                                    sep = "\t",
                                    stringsAsFactors = FALSE,
                                    verbose = FALSE,
                                    data.table = TRUE,
                                    showProgress = TRUE,
                                    header = TRUE,
                                    fill = TRUE,
                                    skip = "Hugo_Symbol",
                                    quote = "")
  }
  ids = coding_ssms %>%
    pull(Tumor_Sample_Barcode) %>%
    unique()

  return(ids)
}
