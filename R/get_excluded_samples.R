#' @title Get Excluded Samples.
#'
#' @description Exclude samples that have been excluded from certain
#' analyses and drop from merges.
#'
#' @details Specify the tool or pipeline responsible for generating the
#' files with
#' `tool_name` and this function will return a vector of excluded sample IDs.
#'
#' @param tool_name The tool or pipeline that generated the files
#' (should be the same for all). Default is "slms-3".
#'
#' @return A vector of sample IDs.
#'
#' @import readr dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#'   excluded_samp <- get_excluded_samples()
#' }
#' 
#' @keywords internal
get_excluded_samples <- function(tool_name = "slms-3") {
  base <- check_config_and_value("repo_base")

  # check for missingness
  path <- paste0(base, "config/exclude.tsv")
  if (!file.exists(path)) {
    print(paste("missing: ", path))
    message("Have you cloned the GAMBL repo and added the path to
    this directory under the local section of your config?")
  }

  excluded_df <- suppressMessages(read_tsv(paste0(base,
                                          "config/exclude.tsv"),
                                  progress = FALSE))
  excluded_samples <- dplyr::filter(excluded_df,
                                    pipeline_exclude == tool_name) %>%
    pull(sample_id)

  return(excluded_samples)
}
