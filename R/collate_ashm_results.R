#' @title Collate ASHM Results.
#'
#' @description Determine the hypermutation status of a few genes.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.results::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return Samples table.
#'
#' @import dplyr tidyr tibble
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#'   sample_table = collate_ashm_results(sample_table = sample_table)
#' }
#' @keywords internal
collate_ashm_results = function(sample_table,
                                seq_type_filter = "genome"){

  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>%
      dplyr::filter(seq_type == seq_type_filter)
      dplyr::select(sample_id, patient_id, biopsy_id)
  }else{
    sample_table = dplyr::filter(sample_table, seq_type == seq_type_filter)
  }
  #just annotate BCL2, MYC and CCND1 hypermutation
  regions_df = data.frame(name = c("CCND1","BCL2","MYC"),
  region = c("chr11:69455000-69459900", "chr18:60983000-60989000", "chr8:128747615-128751834"))
  region_mafs = lapply(regions_df$region, function(x){get_ssm_by_region(region = x, streamlined = FALSE, these_samples_metadata = sample_table)})
  tibbled_data = tibble(region_mafs, region_name = regions_df$name)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)

  unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
    dplyr::select(start, sample_id, region_name)

  tallied = unlisted_df %>%
    group_by(sample_id, region_name) %>%
    tally() %>%
    pivot_wider(values_from = n, names_from = region_name, values_fill = 0, names_prefix = "ashm_")

  sample_table = left_join(sample_table, tallied, by = "sample_id")
}
