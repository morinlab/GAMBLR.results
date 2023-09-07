check_expected_outputs = function(tool_name="battenberg",seq_type_filter="genome"){
  projection = get_template_wildcards("projections")
  #drop irrelevant rows of the metadata based on the scope of the tool etc
  if(tool_name=="battenberg"){
    template_path = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$cnv$battenberg)
    extra_wildcards = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$cnv$battenberg_wildcards)
    #in the current setup, this drops unmatched samples (could be hard-coded but using the config is more transparent)
    relevant_metadata = dplyr::filter(all_meta,base::get(names(extra_wildcards)[1]) == unname(extra_wildcards[1]))
    runs = dplyr::filter(relevant_metadata,seq_type == seq_type_filter) %>%
      mutate(tumour_sample_id = sample_id)
  }else if(tool_name=="slms_3"){
    runs = get_runs_table(seq_type_filter = seq_type_filter) %>% mutate(pair_status = pairing_status)
    vcf_base = get_template_wildcards("vcf_base_name")
    runs = mutate(runs,vcf_base_name = vcf_base)
    #runs = mutate(runs,target_builds = projection)
    runs = expand_grid(runs,target_builds=projection)
    template_path = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$ssm$template$clustered$deblacklisted)
  }
  w = grob_wildcards(template_path)

  #use the unix group and seq_type from the actual metadata and all available projections
  seq_type = seq_type_filter

  runs_files = mutate(runs,outfile=glue::glue(template_path))

}
