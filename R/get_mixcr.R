#' @title Get MiXCR.
#'
#' @description Retrieve the top (most abundant) MiXCR clone for each BCR and/or TCR chain of each sample.
#'
#' @details For every sample in `these_samples_metadata` and every chain of the requested
#' `receptor_type`, the per-sample MiXCR clone table is loaded from the 99-outputs directory
#' and the clone with the highest readFraction is kept. BCR chains (IGH, IGK, IGL) are read
#' from the merged clone tables, which include the glycosylation annotations. TCR chains
#' (TRA, TRB, TRD, TRG) have no glycosylation output, so `num_glycosylation_sites_cdr` and
#' `potential_ags` are NA for those rows. Chains with a missing or empty clone file for a
#' sample are skipped.
#'
#' @param these_samples_metadata A metadata table with (at least) sample_id, seq_type and unix_group columns.
#' @param receptor_type Which receptor chains to retrieve. One of "BCR" (IGH, IGK, IGL), "TCR" (TRA, TRB, TRD, TRG) or "both". Default is "BCR".
#' @param base_path Optional base directory to override project_base from the config.
#' @param verbose Set to TRUE to print the path of each file being read. Default is FALSE.
#'
#' @return A data frame with one row per sample/chain combination with the columns
#' sample_id, seq_type, chain, cloneId, readFraction, bestVHit, bestDHit, bestJHit, bestCHit,
#' num_glycosylation_sites_cdr and potential_ags.
#'
#' @import dplyr readr glue purrr GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' my_meta = get_gambl_metadata(seq_type_filter = "mrna")
#' bcr_clones = get_mixcr(these_samples_metadata = my_meta)
#'
#' #BCR and TCR clones
#' all_clones = get_mixcr(these_samples_metadata = my_meta,
#'                        receptor_type = "both")
#' }
#'
get_mixcr = function(these_samples_metadata,
                     receptor_type = "BCR",
                     base_path,
                     verbose = FALSE){

  if(missing(these_samples_metadata)){
    stop("Please provide a metadata table with these_samples_metadata")
  }
  receptor_type = match.arg(receptor_type, c("BCR", "TCR", "both"))
  bcr_chains = c("IGH", "IGK", "IGL")
  tcr_chains = c("TRA", "TRB", "TRD", "TRG")
  chains = switch(receptor_type,
                  BCR = bcr_chains,
                  TCR = tcr_chains,
                  both = c(bcr_chains, tcr_chains))

  if(missing(base_path)){
    base_path = GAMBLR.helpers::check_config_and_value("project_base")
  }
  merged_template = paste0(base_path,
                           GAMBLR.helpers::check_config_and_value("results_flatfiles$mixcr$merged"))
  clones_template = paste0(base_path,
                           GAMBLR.helpers::check_config_and_value("results_flatfiles$mixcr$clones"))

  sample_groups = these_samples_metadata %>%
    dplyr::select(sample_id, seq_type, unix_group) %>%
    distinct()

  combos = tidyr::expand_grid(sample_groups, chain = chains)

  get_top_clone = function(sample_id, seq_type, unix_group, chain){
    if(chain %in% bcr_chains){
      clone_path = glue::glue(merged_template)
    }else{
      clone_path = glue::glue(clones_template)
    }
    #empty files mean no clones were assembled for that chain
    if(!file.exists(clone_path) || file.size(clone_path) == 0){
      return(NULL)
    }
    if(verbose){
      message(paste("loading from:", clone_path))
    }
    clones = suppressMessages(read_tsv(clone_path,
                                       col_types = cols(.default = "c"),
                                       progress = FALSE))
    if(nrow(clones) == 0){
      return(NULL)
    }
    if(!"num_glycosylation_sites_cdr" %in% colnames(clones)){
      clones = mutate(clones,
                      num_glycosylation_sites_cdr = NA_character_,
                      potential_ags = NA_character_)
    }
    clones %>%
      mutate(sample_id = sample_id,
             seq_type = seq_type,
             chain = chain,
             readFraction = as.numeric(readFraction)) %>%
      slice_max(readFraction, n = 1, with_ties = FALSE) %>%
      dplyr::select(sample_id, seq_type, chain, cloneId, readFraction, starts_with("best"),
                    num_glycosylation_sites_cdr, potential_ags)
  }

  top_clones = purrr::pmap(combos, get_top_clone) %>%
    bind_rows()
  if(nrow(top_clones) == 0){
    message("No MiXCR clones found for the supplied samples")
    return(top_clones)
  }
  top_clones = top_clones %>%
    mutate(cloneId = as.integer(cloneId),
           num_glycosylation_sites_cdr = as.integer(num_glycosylation_sites_cdr)) %>%
    arrange(sample_id, seq_type, chain)
  return(top_clones)
}
