#' @title Create cBioPortal Study.
#'
#' @description Wrapper function for creating a import-ready cBioPortal study.
#'
#' @details This function internally calls [GAMBLR::setup_study], [GAMBLR::setup_fusions], [GAMBLR::finalize_study] and [GAMBLR::study_check] to generate all necessary files for importing a study into cBioPortal.
#' This function was developed to streamline this step and at the same time ensure that the study information and selected data type is consistent throughout the individual steps of generating a study.
#' In addition, the user can also control if the generated study should be checked for sample IDs in case lists that are not described in the clinical file.
#' This potentially will prevent an annoying error that prevents the study to be imported into the active cBioPortal instance, default is TRUE.
#' Fusions are also handled based on the selected seq type (`this_seqtype`).
#'
#' @param this_seqtype The seq type you want to generate a study for. Default is "genome".
#' @param short_name A concise name for your portal project. Default is "GAMBL".
#' @param human_friendly_name A slightly more verbose name for your project. Default is "GAMBL data".
#' @param project_name Unique ID for your project. Default is "gambl_all".
#' @param description A verbose description of your data set. This is what the study will be named when accessing it through cBioPortal. Default is "GAMBL data from genome".
#' @param gambl_maf MAF origin.
#' @param gambl_icgc_maf ICGC MAF origin.
#' @param cancer_type Cancer types included in study, default is "mixed".
#' @param overwrite Flag to specify that files should be overwritten if they exist. Default is TRUE.
#' @param check_study Boolean parameter that controls if the generated study should be checked for sample IDs in case lists, that are not described in the clinical file. Default is TRUE.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return Nothing. Rather, this function generates all files necessary for successfully importing a study into an active cBioPortal instance.
#'
#' @examples
#' \dontrun{
#' #generate cBioPortal study for all GAMBL genome samples:
#' cbioportal_create()
#'
#' #generate a cBioPortal study for all GAMVL capture samples:
#' cbioportal_create(this_seqtype = "capture", description = "GAMBL data from exomes")
#' }
#' 
cbioportal_create = function(this_seqtype = "genome",
                             short_name = "GAMBL",
                             human_friendly_name = "GAMBL data",
                             project_name = "gambl_genome",
                             description = "GAMBL data from genome",
                             gambl_maf = "maf_slms3_hg19",
                             gambl_icgc_maf = "maf_slms3_hg19_icgc",
                             cancer_type = "mixed",
                             overwrite = TRUE,
                             check_study = TRUE,
                             out_dir){

  #setup study
  ids = setup_study(seq_type_filter = this_seqtype,
                    short_name = short_name,
                    human_friendly_name = human_friendly_name,
                    project_name = project_name,
                    description = description,
                    overwrite = overwrite,
                    out_dir = out_dir)

  #setup fusions
  if(this_seqtype == "genome"){
    fusion_ids = setup_fusions(short_name = short_name,
                               human_friendly_name = human_friendly_name,
                               project_name = project_name,
                               description = description,
                               gambl_maf = gambl_maf,
                               gambl_icgc_maf = gambl_icgc_maf,
                               out_dir = out_dir)

  }else if(this_seqtype == "capture"){
    message("The selected seq type is capture, no fusions will be generated...")

  }else{
    stop("Please enter a valid seq_type (i.e genome or capture)...")
  }

  #finalize study
  if(this_seqtype == "genome"){
    finalize_study(seq_type_filter = this_seqtype,
                   short_name = short_name,
                   human_friendly_name = human_friendly_name,
                   project_name = project_name,
                   description = description,
                   cancer_type = cancer_type,
                   these_sample_ids = c(ids, fusion_ids),
                   overwrite = overwrite, out_dir = out_dir)

  }else if(this_seqtype == "capture"){
    finalize_study(seq_type_filter = this_seqtype,
                   short_name = short_name,
                   human_friendly_name = human_friendly_name,
                   project_name = project_name,
                   description = description,
                   cancer_type = cancer_type,
                   these_sample_ids = ids,
                   overwrite = overwrite, out_dir = out_dir)

  }else{
    stop("Please enter a valid seq_type (i.e genome or capture)...")
  }

  #check study
  study_check(data_clinical_samples_path = "data_clinical_samples.txt",
              data_fusions_path = "data_fusions.txt",
              cases_fusions_path = "case_lists/cases_fusion.txt",
              cases_all_path = "case_lists/cases_all.txt",
              cases_sequenced_path = "case_lists/cases_sequenced.txt",
              project_name = project_name,
              out_dir = out_dir)
}
