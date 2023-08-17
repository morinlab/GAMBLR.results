#' @title Collate PGA results for samples with CN data.
#'
#' @description Expand a metadata table horizontally with PGA metrics.
#'
#' @details Helper function called by `collate_results`, not meant for out-of-package usage.
#'
#' @param these_samples_metadata The metadata to be expanded with sample_id column.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @noRd
#'
#' @return data frame
#' @import dplyr
#'
#' @examples
#' # For genomes
#' meta <- get_gambl_metadata()
#' pga_metrics <- collate_pga(these_samples_metadata = meta)
#' # For exomes
#' meta_capture <- get_gambl_metadata(seq_type_filter = "capture")
#' pga_metrics_capture <- collate_pga(these_samples_metadata = meta_capture)
#'
collate_pga <- function(
    these_samples_metadata,
    this_seq_type = "genome"
) {

    message(
        "Collating the PGA results ..."
    )
    # Currently only works for genomes
    if(! this_seq_type %in% c("genome", "capture")) {
        stop("Please provide a valid seq_type (\"genome\" or \"capture\").")
    }

    # Default to all samples if sample table is missing
    if (missing(these_samples_metadata)) {
        message("No sample table was provided. Defaulting to all metadata ...")
        these_samples_metadata <- get_gambl_metadata(
            seq_type_filter = this_seq_type
        )
    }

    # Get the CN segments
    multi_sample_seg <- get_sample_cn_segments(
        sample_list = these_samples_metadata$sample_id,
        multiple_samples = TRUE,
        this_seq_type = this_seq_type
    ) %>%
    dplyr::rename("sample" = "ID")

    these_samples_pga <- calculate_pga(
        this_seg = multi_sample_seg
    )

    these_samples_metadata <- left_join(
        these_samples_metadata,
        these_samples_pga
    )

    return(these_samples_metadata)

}
