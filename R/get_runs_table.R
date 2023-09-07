get_runs_table = function(seq_type_filter="genome"){
  t_meta = get_gambl_metadata(tissue_status_filter = c("tumour"),seq_type_filter=seq_type_filter) %>%
    dplyr::select(sample_id,patient_id,seq_type,genome_build,pairing_status,unix_group) %>%
    rename("tumour_sample_id"="sample_id")
  n_meta = get_gambl_metadata(tissue_status_filter = c("normal"),seq_type_filter=seq_type_filter) %>%
    dplyr::select(sample_id,patient_id,seq_type,genome_build) %>% rename("normal_sample_id"="sample_id")
  runs_df = left_join(t_meta,n_meta,by=c("patient_id","seq_type","genome_build"))
  #fill in normal_sample_id for unmatched cases
  unmatched_df = get_unmatched_normals(seq_type_filter=seq_type_filter)
  runs_df = left_join(runs_df,unmatched_df,by=c("seq_type","genome_build","unix_group")) %>%
    mutate(normal_sample_id=ifelse(is.na(normal_sample_id.x),normal_sample_id.y,normal_sample_id.x)) %>%
    select(-normal_sample_id.x, -normal_sample_id.y)
  return(runs_df)
}
