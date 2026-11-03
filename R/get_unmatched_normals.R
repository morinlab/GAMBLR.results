get_unmatched_normals = function(seq_type_filter){
  a = GAMBLR.helpers::check_config_value(config::get("unmatched_normal_ids"))
  #flatten unix_group -> seq_type -> genome_build into one row per normal
  df = tibble::tibble(unix_group = names(a), normal_sample_id = a) %>%
    tidyr::unnest_longer(normal_sample_id, indices_to = "seq_type") %>%
    tidyr::unnest_longer(normal_sample_id, indices_to = "genome_build") %>%
    dplyr::mutate(normal_sample_id = as.character(normal_sample_id)) %>%
    dplyr::filter(seq_type == seq_type_filter)
  return(df)
}
