#' @title Compare Mutation Flavour.
#'
#' @description Get a MAF that is just the variants unique to one of two flavours of variant calls available.
#'
#' @details Subset a MAF to only have variants that are unique to one flavour (specified with `flavour1`).
#' This function is currently not exported, since there is only one flavour available at the moment (see docs for [GAMBLR::get_ssm_by_sample]).
#'
#' @param these_sample_ids A vector of sample IDs to be included.
#' @param flavour1 First flavour of variant calls, to be returned as unique if not present in flavour2. Default is "clustered".
#' @param flavour2 Second flavour of variant calls.
#'
#' @return a list with MAFs that are only present in flavour1.
#'
#' @noRd
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
compare_mutation_flavour = function(these_sample_ids,
                                    flavour1 = "clustered",
                                    flavour2 = ""){

  these_dfs = list()
  for(this_sample_id in these_sample_ids){
    message(this_sample_id)
    maf1 = get_ssm_by_sample(this_sample_id, flavour = flavour1, this_seq_type = seq_type)
    maf2  = get_ssm_by_sample(this_sample_id, flavour = flavour2)
    maf1_only = intersect_maf(maf1, maf2)
    these_dfs[[this_sample_id]] = maf1_only
  }
  this_maf = rbindlist(these_dfs, use.names = TRUE)
  return(this_maf)
}
