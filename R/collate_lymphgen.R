#' @title Compute LymphGen results for a sample scope.
#'
#' @description Core computation behind `collate_lymphgen()`, extracted so
#' it can also be called directly by [collate_results_db()] for just the
#' subset of samples missing from its cache table. Note that unlike
#' `compute_ssm_results_core()`, the underlying LymphGen flavour files get
#' loaded in full regardless of scope -- they're small enough that this
#' isn't worth avoiding -- `these_samples_metadata`, when provided, only
#' narrows what's returned, not what's loaded.
#'
#' @param these_samples_metadata Optional parameter with metadata filtered
#' for sample_ids of interest. If provided, the result is filtered down to
#' just those sample_ids.
#' @param lymphgen_version Version of selected lymphgen, default is "default".
#' @param tidy Boolean parameter, set to TRUE for tidy format (i.e long
#' format with no columns dropped). Default is FALSE, which returns the
#' data in a wide format, keeping both the original Subtype. Prediction and
#' tidied LymphGen values and puts the values from each "flavour" in its
#' own column.
#'
#' @return A df with lymphgen information, keyed by `sample_id`.
#'
#' @import dplyr tidyr readr stringr glue GAMBLR.helpers
#'
#' @keywords internal
#' @noRd
compute_lymphgen_core <- function(these_samples_metadata,
                                  lymphgen_version = "default",
                                  tidy = FALSE){

  #TODO Update the key in the config to match the version once updated, as discussed on PR.
  if(lymphgen_version == "default"){
    lymphgen_template = GAMBLR.helpers::check_config_value(config::get("results_versioned")$lymphgen_template$default)
  }else{
    stop("Currently, only lymphgen_version = default is accepted")
  }

  #repo base
  repo_base = GAMBLR.helpers::check_config_value(config::get("repo_base"))
  flavours = GAMBLR.helpers::check_config_value(config::get("results_merged_wildcards")$lymphgen_template)
  flavour = str_split(flavours, pattern = ",")
  flavour = unlist(flavour)
  lymphgen_path = paste0(repo_base, lymphgen_template)

  load_lymphgen = function(flavour, lymphgen_path){
    lg_path = glue::glue(lymphgen_path)
    if(!file.exists(lg_path)){ #ignore missing flavours.
      return()
    }
    lg_df = suppressMessages(read_tsv(lg_path)) %>%
      mutate(flavour = flavour) #append the flavour in its own column called "flavour".
    return(lg_df)
  }

  lymphgen_results = lapply(flavour, load_lymphgen, lymphgen_path = lymphgen_path)
  lymphgen_results = bind_rows(lymphgen_results) #get lymphgen results tables stacked on top of each other, with the results from each flavour identified by the `flavour` column.
  lymphgen_results = tidy_lymphgen(lymphgen_results, lymphgen_column_in = "Subtype.Prediction", lymphgen_column_out = "LymphGen")
  colnames(lymphgen_results)[1] = "sample_id"

  if(!tidy){
    result = lymphgen_results %>%
      select(sample_id, Subtype.Prediction, LymphGen, flavour) %>%
      pivot_wider(names_from = flavour,
                  values_from = c(Subtype.Prediction, LymphGen),
                  names_glue = "{.value}_{flavour}")
  }else{
    result = lymphgen_results
  }

  if(!missing(these_samples_metadata)){
    # every requested sample_id gets a row back, even one absent from every
    # loaded flavour file (NA in the derived columns), not just those that
    # happened to have a classification. Without this, collate_results_db()
    # would never see such a sample as "already computed" -- its cache table
    # would never gain a row for it -- so every future run would reload and
    # re-scan all flavour files for it again, forever, instead of caching
    # the "no classification available" result the way a real NA does.
    result = dplyr::select(these_samples_metadata, sample_id) %>%
      dplyr::left_join(result, by = "sample_id")
  }

  return(result)
}

#' @title Collate Lymphgen.
#'
#' @description Expand a sample_table (metadata) horizontally with different flavours of lymphgen data.
#'
#' @details This function takes a sample table (metadata) and adds different flavours of lymphgen data.
#' It is possible to call this function with an already subset metadata table (with sample IDs of interest) with `these_samples_metadata`.
#' If this is done, the function will join the lymphgen data with this table. Currently, the only supported `lymphgen_version` is "default".
#' For more information refer to the function examples.
#'
#' @param these_samples_metadata Optional parameter with metadata filtered for sample_ids of interest. If provided, this function will join lymphgen with this metadata, regardless of tidy TRUE/FALSE.
#' @param lymphgen_version Version of selected lymphgen, default is "default".
#' @param tidy Boolean parameter, set to TRUE for tidy format (i.e long format with no columns dropped). Default is FALSE, which returns the data in a wide format, keeping both the original Subtype. Prediction and tidied LymphGen values and puts the values from each "flavour" in its own column.
#'
#' @return A df with lymphgen information.
#'
#' @import dplyr tidyr readr stringr glue GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' this_meta = get_gambl_metadata()
#' dlbcl_meta = dplyr::filter(this_meta, pathology == "DLBCL")
#'
#' wide_lymphgen = collate_lymphgen(these_samples_metadata = dlbcl_meta,
#'                                  lymphgen_version = "default",
#'                                  tidy = FALSE)
#'}
collate_lymphgen = function(these_samples_metadata,
                            lymphgen_version = "default",
                            tidy = FALSE){

  if(missing(these_samples_metadata)){
    result = compute_lymphgen_core(lymphgen_version = lymphgen_version, tidy = tidy)
    return(result)
  }

  new_cols = compute_lymphgen_core(these_samples_metadata = these_samples_metadata,
                                   lymphgen_version = lymphgen_version,
                                   tidy = tidy)
  result = left_join(these_samples_metadata, new_cols, by = "sample_id")

  return(result)
}
