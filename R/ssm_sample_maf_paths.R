# Paths of one sample's SSM MAFs, as get_ssm_by_sample() reads them: the deblacklisted and augmented MAFs, on the
# GSC (config "default") and in the local mirror (the active config's project_base). Shared by get_ssm_by_sample()
# and the collate cache's provenance (ssm_input_files()), so both always point at the same files.
#
# these_samples_metadata: one row (sample_id, seq_type, unix_group, genome_build, pairing_status, normal_sample_id)
# Returns a list with full_maf_path, local_full_maf_path, aug_maf_path, local_aug_maf_path and vcf_base_name.
.ssm_sample_maf_paths = function(these_samples_metadata, projection = "grch37", flavour = "clustered"){
  if(!flavour %in% c("clustered", "sage")){
    stop("flavour must be \"clustered\" or \"sage\"")
  }
  # variables the config path templates refer to (glue resolves them in this frame)
  seq_type = dplyr::pull(these_samples_metadata, seq_type)
  sample_id = dplyr::pull(these_samples_metadata, sample_id)
  tumour_sample_id = sample_id
  unix_group = dplyr::pull(these_samples_metadata, unix_group)
  genome_build = dplyr::pull(these_samples_metadata, genome_build)
  target_builds = projection
  pair_status = dplyr::pull(these_samples_metadata, pairing_status)
  # Unmatched samples were called against the unmatched normal in the GAMBLR config; keep that up to date if/when
  # any new references are added.
  if(pair_status == "unmatched"){
    keys = paste0("unmatched_normal_ids$", unix_group, "$", seq_type, "$", genome_build)
    normal_sample_id = check_config_and_value(keys)
  }else{
    normal_sample_id = dplyr::pull(these_samples_metadata, normal_sample_id)
  }
  base_path = ""
  vcf_base_name = if(flavour == "clustered") "slms-3.final" else "sage.combined"

  remote_base = check_config_and_value("project_base", config_name = "default")
  local_base = check_config_and_value("project_base")
  template_path = function(which){
    template = check_config_and_value(paste0("results_flatfiles$ssm$template$", flavour, "$", which),
                                      config_name = "default")
    unname(unlist(glue::glue(template)))
  }
  deblacklisted = template_path("deblacklisted")
  augmented = template_path("augmented")
  list(full_maf_path = paste0(remote_base, deblacklisted),
       local_full_maf_path = paste0(local_base, deblacklisted),
       aug_maf_path = paste0(remote_base, augmented),
       local_aug_maf_path = paste0(local_base, augmented),
       vcf_base_name = vcf_base_name)
}
