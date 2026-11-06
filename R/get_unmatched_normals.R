get_unmatched_normals = function(seq_type_filter, required = NULL){
  a = GAMBLR.helpers::check_config_value(config::get("unmatched_normal_ids"))
  #flatten unix_group -> seq_type -> genome_build into one row per normal
  df = tibble::tibble(unix_group = names(a), normal_sample_id = a) %>%
    tidyr::unnest_longer(normal_sample_id, indices_to = "seq_type") %>%
    tidyr::unnest_longer(normal_sample_id, indices_to = "genome_build") %>%
    dplyr::mutate(normal_sample_id = as.character(normal_sample_id)) %>%
    dplyr::filter(seq_type == seq_type_filter)
  #required: samples (unix_group, seq_type, genome_build) that need an unmatched normal
  if(!is.null(required)){
    missing_combos = required %>%
      dplyr::count(unix_group, seq_type, genome_build, name = "n_samples") %>%
      dplyr::anti_join(df, by = c("unix_group", "seq_type", "genome_build"))
    if(nrow(missing_combos) > 0){
      stop(
        "No unmatched_normal_ids entry configured for these unix_group/seq_type/genome_build combinations:\n",
        paste0("  ", missing_combos$unix_group, "/", missing_combos$seq_type, "/", missing_combos$genome_build,
               " (", missing_combos$n_samples, " samples)", collapse = "\n"),
        "\nAdd them to the unmatched_normal_ids config, or check these samples' metadata for errors."
      )
    }
  }
  return(df)
}
