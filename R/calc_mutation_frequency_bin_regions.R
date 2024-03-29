#' @title Mutation counts across sliding windows for multiple regions.
#'
#' @description Obtain a long tidy or wide matrix of mutation counts across sliding windows for multiple regions.
#'
#' @details This function takes a metadata table with `these_samples_metadata` parameter and internally calls [GAMBLR.results::calc_mutation_frequency_bin_region] (that internally calls [GAMBLR.results::get_ssm_by_regions]).
#' to retrieve mutation counts for sliding windows across one or more regions. May optionally provide any combination of a maf data frame, existing metadata, or a regions data frame or named vector.
#' The heatmap plotting portion of this function has been moved to [GAMBLR.results::heatmap_mutation_frequency_bin].
#'
#' @param regions_list Named vector of regions in the format c(name1 = "chr:start-end", name2 = "chr:start-end"). If neither regions nor regions_bed is specified, the function will use GAMBLR aSHM region information.
#' @param regions_bed Data frame of regions with four columns (chrom, start, end, name).
#' @param these_samples_metadata Metadata with at least sample_id column. If not providing a maf data frame, seq_type is also required.
#' @param these_sample_ids Vector of sample IDs. Metadata will be subset to sample IDs present in this vector.
#' @param this_seq_type Optional vector of seq_types to include in heatmap. Default c("genome", "capture"). Uses default seq_type priority for samples with >1 seq_type.
#' @param maf_data Optional maf data frame. Will be subset to rows where Tumor_Sample_Barcode matches provided sample IDs or metadata table. If not provided, maf data will be obtained with get_ssm_by_regions().
#' @param region_padding Amount to pad the start and end coordinates by. Default 1000.
#' @param projection Genome build the function will operate in. Ensure this matches your provided regions and maf data for correct chr prefix handling. Default "grch37".
#' @param drop_unmutated Whether to drop bins with 0 mutations. If returning a matrix format, this will only drop bins with no mutations in any samples.
#' @param skip_regions Optional character vector of genes to exclude from the default aSHM regions.
#' @param only_regions Optional character vector of genes to include from the default aSHM regions.
#' @param slide_by Slide size for sliding window. Default 100.
#' @param window_size Size of sliding window. Default 500.
#' @param return_format Return format of mutations. Accepted inputs are "long" and "wide". Long returns a data frame of one sample ID/window per row. Wide returns a matrix with one sample ID per row and one window per column. Using the "wide" format will retain all samples and windows regardless of the drop_unmutated or min_count_per_bin parameters. Default wide.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flat files (only works for streamlined data, not full MAF details).
#' @param mode Only works with indexed flat files. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#' @return A table of mutation counts for sliding windows across one or more regions. May be long or wide.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr tibble
#' @export
#'
#' @examples
#' #load metadata.
#' metadata = get_gambl_metadata()
#' dlbcl_bl_meta = dplyr::filter(metadata, pathology %in% c("DLBCL", "BL"))
#'
#' #bring together all derived sample-level results from many GAMBL pipelines.
#' dlbcl_bl_meta = collate_results(join_with_full_metadata = TRUE,
#'                                 these_samples_metadata = dlbcl_bl_meta)
#'
#' #get ashm regions
#' some_regions = GAMBLR.data::grch37_ashm_regions
#'
#' mut_count_matrix <- calc_mutation_frequency_bin_region(
#'    these_samples_metadata = dlbcl_bl_meta,
#'    regions_bed = some_regions
#' )
#'
calc_mutation_frequency_bin_regions <- function(
  regions_list = NULL,
  regions_bed = NULL,
  these_samples_metadata = NULL,
  these_sample_ids = NULL,
  this_seq_type = c("genome", "capture"),
  maf_data = NULL,
  projection = "grch37",
  region_padding = 1000,
  drop_unmutated = FALSE,
  skip_regions = NULL,
  only_regions = NULL,
  slide_by = 100,
  window_size = 500,
  return_format = "wide",
  from_indexed_flatfile = TRUE,
  mode = "slms-3"
) {
  regions <- process_regions(
    regions_list = regions_list,
    regions_bed = regions_bed,
    region_padding = region_padding,
    skip_regions = skip_regions,
    only_regions = only_regions
  )
  regions_bed <- regions$regions_bed
  regions <- regions$regions_list

  if (
    (str_detect(regions_bed$chrom[1], "chr") & projection == "grch37") |
      (!str_detect(regions_bed$chrom[1], "chr") & projection == "hg38")
  ) {
    stop("chr prefixing status of provided regions and specified projection don't match. ")
  }
  # Harmonize metadata and sample IDs
  get_meta <- id_ease(
    these_samples_metadata,
    these_sample_ids,
    this_seq_type
  )
  metadata <- get_meta
  these_sample_ids <- get_meta$sample_id %>% unique

  # Obtain sliding window mutation frequencies for all regions
  dfs <- parallel::mclapply(names(regions), function(x) {
    df <- calc_mutation_frequency_bin_region(
      region = regions[x],
      these_samples_metadata = metadata,
      maf_data = maf_data,
      projection = projection,
      drop_unmutated = drop_unmutated,
      slide_by = slide_by,
      window_size = window_size,
      min_count_per_bin = 0,
      return_count = TRUE,
      from_indexed_flatfile = from_indexed_flatfile,
      mode = mode
    ) %>%
      dplyr::mutate(name = x)
    return(df)
  })


  all <- dplyr::bind_rows(dfs) %>%
    dplyr::distinct(bin, sample_id, .keep_all = TRUE)

  # If none of the samples are mutated, return the mutation frequency df and exit.
  if (max(all$mutation_count) == 0) {
    message("No mutations found in specified regions for specified samples. Exiting. ")
    return(all)
  }

  if (return_format == "wide") {
    # Convert mutation frequency table to a matrix
    all_wide <- all %>%
      dplyr::select(sample_id, mutation_count, bin) %>%
      pivot_wider(
        names_from = bin,
        values_from = mutation_count,
        values_fill = 0
      )
    return(all_wide)
  } else {
    return(all)
  }
}
