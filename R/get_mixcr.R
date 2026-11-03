#' @title Get MiXCR.
#'
#' @description Retrieve the top (most abundant) MiXCR clone for each BCR and/or TCR chain of each sample.
#'
#' @details For every sample in `these_samples_metadata` and every requested chain
#' (from `chains`, or else all chains of `receptor_type`), the per-sample dominant MiXCR clone table (with glycosylation annotations)
#' is loaded from the 99-outputs directory and the clone with the highest readFraction is kept.
#' Chains with a missing or empty clone file for a sample are skipped. TCR outputs are only
#' generated when requested, so most samples will not have TCR clones.
#'
#' @param these_samples_metadata A metadata table with (at least) sample_id, seq_type and unix_group columns.
#' @param receptor_type Which receptor chains to retrieve. One of "BCR" (IGH, IGK, IGL), "TCR" (TRA, TRB, TRD, TRG) or "both". Default is "BCR". Ignored if `chains` is provided.
#' @param chains Optional vector of specific chains to retrieve, e.g. "IGH" or c("IGH", "IGK"). Overrides `receptor_type`.
#' @param base_path Optional base directory to override project_base from the config.
#' @param num_cores Number of cores used to read the clone files in parallel with [parallel::mclapply]. Default is 8. Set to 1 to read sequentially.
#' @param verbose Set to TRUE to print the path of each file being read. Default is FALSE.
#'
#' @return A data frame with one row per sample/chain combination with the columns
#' sample_id, seq_type, chain, cloneId, readFraction, bestVHit, bestDHit, bestJHit, bestCHit,
#' num_glycosylation_sites_cdr and potential_ags.
#'
#' @import dplyr readr glue parallel GAMBLR.helpers
#' @export
#'
#' @examples
#' 
#' my_meta = get_gambl_metadata() %>% 
#'  dplyr::filter(seq_type == "mrna", cohort == "DLBCL_DLC")
#' bcr_clones = get_mixcr(these_samples_metadata = my_meta)
#'
#' #BCR and TCR clones
#' all_clones = get_mixcr(these_samples_metadata = my_meta,
#'                        receptor_type = "both")
#' head(all_clones)
#'
#' #only heavy and kappa chain clones
#' igh_igk_clones = get_mixcr(these_samples_metadata = my_meta,
#'                            chains = c("IGH", "IGK"))
#' head(igh_igk_clones)
#'
get_mixcr = function(these_samples_metadata,
                     receptor_type = "BCR",
                     chains,
                     base_path,
                     num_cores = 8,
                     verbose = FALSE){

  if(missing(these_samples_metadata)){
    stop("Please provide a metadata table with these_samples_metadata")
  }
  bcr_chains = c("IGH", "IGK", "IGL")
  tcr_chains = c("TRA", "TRB", "TRD", "TRG")
  if(missing(chains)){
    receptor_type = match.arg(receptor_type, c("BCR", "TCR", "both"))
    chains = switch(receptor_type,
                    BCR = bcr_chains,
                    TCR = tcr_chains,
                    both = c(bcr_chains, tcr_chains))
  }else{
    chains = unique(toupper(chains))
    invalid_chains = setdiff(chains, c(bcr_chains, tcr_chains))
    if(length(invalid_chains) > 0){
      stop(paste0("Invalid chain(s): ", paste(invalid_chains, collapse = ", "),
                  ". Valid chains are: ", paste(c(bcr_chains, tcr_chains), collapse = ", ")))
    }
  }

  if(missing(base_path)){
    base_path = GAMBLR.helpers::check_config_and_value("project_base")
  }
  clone_template = paste0(base_path,
                          GAMBLR.helpers::check_config_and_value("results_flatfiles$mixcr"))

  sample_groups = these_samples_metadata %>%
    dplyr::select(sample_id, seq_type, unix_group) %>%
    distinct()

  combos = tidyr::expand_grid(sample_groups, chain = chains) %>%
    mutate(clone_path = as.character(glue::glue_data(., clone_template)))

  get_top_clone = function(sample_id, seq_type, unix_group, chain, clone_path){
    #checked here so the slow filesystem calls run in parallel.
    #file.size is NA for missing files; empty files mean no clones were assembled for that chain
    clone_size = file.size(clone_path)
    if(is.na(clone_size) || clone_size == 0){
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
    clones %>%
      mutate(sample_id = sample_id,
             seq_type = seq_type,
             chain = chain,
             readFraction = as.numeric(readFraction)) %>%
      slice_max(readFraction, n = 1, with_ties = FALSE) %>%
      dplyr::select(sample_id, seq_type, chain, cloneId, readFraction, starts_with("best"),
                    num_glycosylation_sites_cdr, potential_ags)
  }

  top_clones = parallel::mclapply(seq_len(nrow(combos)), function(i){
    do.call(get_top_clone, as.list(combos[i, ]))
  }, mc.cores = num_cores)

  failed = vapply(top_clones, function(x){inherits(x, "try-error")}, logical(1))
  if(any(failed)){
    stop(paste("Failed to read:", paste(combos$clone_path[failed], collapse = ", "),
               "\n", top_clones[failed][[1]]))
  }
  top_clones = bind_rows(top_clones)
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
