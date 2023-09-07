add_prps_result = function(incoming_metadata){
  prps_res = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/derived_and_curated_metadata/outputs/BL_dhitsig_PRPS.tsv"))
  colnames(prps_res)[1] = "sample_id"
  prps_res = dplyr::select(prps_res, sample_id, PRPS_score, PRPS_class)

  #need to associate each sample with a patient ID then annotate the metadata based on patient ID
  patient_meta_g = get_gambl_metadata(seq_type_filter = "genome") %>%
    dplyr::select(sample_id, patient_id)

  patient_meta_r = get_gambl_metadata(seq_type_filter = "mrna") %>%
    dplyr::select(sample_id, patient_id)

  patient_meta = bind_rows(patient_meta_g, patient_meta_r)
}
