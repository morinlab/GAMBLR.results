#' @title Get ASHM Count Matrix.
#'
#' @description Prepare a matrix with one row per sample and one column per region using a set of hypermutated regions.
#'
#' @details Values are the number of mutations in that patient in the region.
#'
#' @param regions_bed A bed file with one row for each region.
#' @param maf_data Optionally provide a data frame in the MAF format, otherwise the database will be used.
#' @param these_samples_metadata This is used to complete your matrix. All GAMBL samples will be used by default. Provide a data frame with at least sample_id for all samples if you are using non-GAMBL data.
#' @param seq_type The seq type to return results for.
#' @param from_indexed_flatfile Boolean parameter set to TRUE per default.
#'
#' @return A matrix.
#'
#' @import dplyr tibble
#' @export
#'
#' @examples
#' regions_bed = dplyr::mutate(GAMBLR.data::somatic_hypermutation_locations_GRCh37_v_latest, name = paste(gene, region, sep = "_"))
#'
#' matrix = get_ashm_count_matrix(regions_bed = regions_bed,
#'                                seq_type = "genome")
#'
get_ashm_count_matrix = function(regions_bed,
                                 maf_data,
                                 these_samples_metadata,
                                 seq_type,
                                 from_indexed_flatfile = TRUE){
  if(missing(seq_type)){
    if(missing(these_samples_metadata)){
      stop("Must supply either the seq_type or a metadata data frame from which it can be retrieved")
    }
    seq_type = head(these_samples_metadata) %>% pull(seq_type)
  }
  if(missing(regions_bed)){
    regions_bed = grch37_ashm_regions
  }
  ashm_maf = get_ssm_by_regions(regions_bed = regions_bed,
                                streamlined = TRUE,
                                seq_type=seq_type,
                                maf_data = maf_data,
                                use_name_column = TRUE,
                                from_indexed_flatfile = from_indexed_flatfile)

  ashm_counted = ashm_maf %>%
    group_by(sample_id, region_name) %>%
    tally()

  if(missing(these_samples_metadata)){
    all_meta = get_gambl_metadata(seq_type_filter=seq_type) %>%
      dplyr::select(sample_id)
  }else{
    all_meta = these_samples_metadata %>%
      dplyr::select(sample_id)
  }
  #fill out all combinations so we can get the cases with zero mutations
  eg = expand_grid(sample_id = pull(all_meta, sample_id), region_name = unique(ashm_counted$region_name))
  all_counts = left_join(eg, ashm_counted) %>%
    mutate(n = replace_na(n, 0)) %>%
    unique() #not sure where the duplicates are coming from but its annoying

  all_counts_wide = pivot_wider(all_counts, id_cols = sample_id, names_from = region_name, values_from = n) %>%
    column_to_rownames(var = "sample_id")

  return(all_counts_wide)
}
