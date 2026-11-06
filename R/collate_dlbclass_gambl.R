#' @title Compute GAMBL DLBclass calls for a sample scope.
#'
#' @description Core computation behind `collate_dlbclass_gambl()`, split out
#' so it can also be called directly by [collate_results_db()] for just the
#' subset of samples missing from its cache table. Returns only the new
#' columns, keyed by `sample_id`. The DLBclass results files are read in full;
#' `these_samples_metadata`, when provided, only narrows what is returned.
#'
#' @param these_samples_metadata Optional metadata with the sample_ids of
#' interest. If provided, every requested sample_id gets exactly one row back
#' (NA where DLBclass has no call). If missing, every sample in the results
#' files is returned.
#' @param flavour Which DLBclass runs to collate: "all" (default), or any of
#' "full" (with SVs and CNVs), "no_sv_no_cnv", "no_sv" and "no_cnv". All
#' columns get the flavour as a suffix (e.g. dlbclass_full, dlbclass_confidence_no_cnv).
#' @param return_probabilities Set to TRUE to keep the per-cluster
#' probabilities (C1-C5). Default is FALSE.
#'
#' @return A data frame with `sample_id` and the DLBclass columns.
#'
#' @import dplyr GAMBLR.helpers
#'
#' @keywords internal
#' @noRd
compute_dlbclass_gambl_core <- function(these_samples_metadata,
                                        flavour = "all",
                                        return_probabilities = FALSE){

    # Define flavour names so paths and objects can be operated on correctly
    all_flavours <- c(
        full = "with_sv--with_cnv",
        no_sv_no_cnv = "no_sv--no_cnv",
        no_sv = "no_sv--with_cnv",
        no_cnv = "with_sv--no_cnv"
    )

    if (identical(flavour, "all")) flavour <- names(all_flavours)

    unknown <- setdiff(flavour, names(all_flavours))
    if (length(unknown) > 0) {
        stop("Unknown flavour(s): ", paste(unknown, collapse = ", "),
             ". Use \"all\" or any of: ", paste(names(all_flavours), collapse = ", "))
    }

    # DLBclass results directory
    output_base <- GAMBLR.helpers::check_config_value(config::get("project_base"))
    dlbclass_dir <- file.path(output_base, "all_the_things/dlbclass-1.0/02-dlbclass/gambl-all--2025-10")
    # TODO: add DLBclass results path to config.yml to move away from hardcoded results filepaths here

    # Read each flavour, rename, and combine
    calls_combined <- NULL

    for (this_flavour in flavour){
        file <- file.path(dlbclass_dir, paste0("gambl-all--", all_flavours[[this_flavour]], "_classified_samples.tsv"))
        suffix <- paste0("_", this_flavour) # e.g. dlbclass_full, dlbclass_confidence_no_cnv

        calls <- data.table::fread(file) %>%
            dplyr::as_tibble() %>%
            dplyr::rename(dlbclass = PredictedCluster, sample_id = sample, dlbclass_confidence = Confidence)

        # Remove probabilities if return_probabilities = FALSE
        if (isFALSE(return_probabilities)){
            calls <- dplyr::select(calls, -dplyr::matches("^C[1-5]$"))
        }

        calls <- dplyr::rename_with(calls, ~ paste0(.x, suffix), -sample_id)

        # full_join keeps samples that are missing from any one flavour
        if (is.null(calls_combined)){
            calls_combined <- calls
        } else {
            calls_combined <- dplyr::full_join(calls_combined, calls, by = "sample_id")
        }
    }

    if (anyDuplicated(calls_combined$sample_id)) {
        warning("Some samples are listed more than once in the DLBclass results files, so they will have more than one row")
    }

    # One row per requested sample_id, including samples with no DLBclass call (NA), so collate_results_db() can cache "no call" instead of re-reading the files for them on every run
    if (!missing(these_samples_metadata)){
        calls_combined <- these_samples_metadata %>%
            dplyr::distinct(sample_id) %>%
            dplyr::left_join(calls_combined, by = "sample_id")
    }

    return(calls_combined)
}

#' @title Collate GAMBL DLBclass calls.
#'
#' @description Expand a metadata table horizontally with DLBclass calls from
#' the GAMBL DLBclass runs (one set of columns per flavour).
#'
#' @param these_samples_metadata Optional metadata with the sample_ids of
#' interest. If provided, the DLBclass columns are joined onto it (any
#' existing copies of those columns are replaced). If missing, the calls for
#' every sample in the results files are returned.
#' @param flavour Which DLBclass runs to collate: "all" (default), or any of
#' "full" (with SVs and CNVs), "no_sv_no_cnv", "no_sv" and "no_cnv". All
#' columns get the flavour as a suffix (e.g. dlbclass_full, dlbclass_confidence_no_cnv).
#' @param return_probabilities Set to TRUE to keep the per-cluster
#' probabilities (C1-C5). Default is FALSE.
#'
#' @return The metadata with the DLBclass columns added, or the calls alone
#' if no metadata is given.
#'
#' @import dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' dlbcl_meta <- suppressMessages(get_gambl_metadata()) %>%
#'     dplyr::filter(pathology == "DLBCL")
#'
#' # All four flavours
#' dlbcl_collated <- collate_dlbclass_gambl(dlbcl_meta)
#' dplyr::count(dlbcl_collated, dlbclass_full)
#'
#' # Full flavour only
#' collate_dlbclass_gambl(dlbcl_meta, flavour = "full") %>%
#'     dplyr::select(sample_id, dlbclass_full, dlbclass_confidence_full)
#'
collate_dlbclass_gambl <- function(these_samples_metadata,
                                   flavour = "all",
                                   return_probabilities = FALSE){

    if (missing(these_samples_metadata)){
        return(compute_dlbclass_gambl_core(flavour = flavour,
                                           return_probabilities = return_probabilities))
    }

    new_cols <- compute_dlbclass_gambl_core(these_samples_metadata = these_samples_metadata,
                                            flavour = flavour,
                                            return_probabilities = return_probabilities)

    # Drop any existing copies of these columns so the join doesn't create .x/.y columns
    these_samples_metadata <- these_samples_metadata %>%
        dplyr::select(-dplyr::any_of(setdiff(colnames(new_cols), "sample_id"))) %>%
        dplyr::left_join(new_cols, by = "sample_id")

    return(these_samples_metadata)
}