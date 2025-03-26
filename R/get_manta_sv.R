#' @title Get Manta SVs
#'
#' @description Retrieve Manta SVs for one or many samples
#'
#' @details Retrieve Manta SVs with additional VCF information to allow for
#' filtering of high-confidence variants.
#' To get SV calls for multiple samples, supply a metadata table via
#' `these_samples_metadata` that has been subset to only those samples.
#' The results will be restricted to the sample_ids within that data frame.
#' This function relies on a set of specific internal functions
#' [GAMBLR.results::get_manta_sv_by_samples] (if `from_cache = FALSE`).
#' This function can also restrict the returned breakpoints within a genomic
#' region specified via `region` (in chr:start-end format).
#' Useful filtering parameters are also available, use `min_vaf` to set the
#' minimum tumour VAF for a SV to be returned and `min_score`
#' to set the lowest Manta somatic score for a SV to be returned.
#' In addition, the user can chose to return all variants, even
#' the ones not passing the filter criteria. To do so,
#' set `pass_filters = FALSE` (defaults to TRUE).
#' 
#' ## **Advanced settings (probably not for you)**
#' 
#' Is it advised to leave the default `from_cache` setting to TRUE.
#' To ensure manta results arre pulled from a pre-generated merge
#' (i.e. the cached result).
#' If set to FALSE in combination with `write_to_file = TRUE`,
#' the function will (re)generate new merged manta calls, if the user has
#' the required file permissions.
#' Note, that if `write_to_file` is set to TRUE, the function defaults
#' `from_cache = FALSE` to avoid nonsense parameter combinations.
#' Is this function not what you are looking for? You may want:
#' [GAMBLR.results::get_combined_sv]
#' After running this or [get_combined_sv], you most likely want to
#' annotate the result using [GAMBLR.utils::annotate_sv]
#'
#' @param these_samples_metadata A metadata data frame to limit the
#' result to sample_ids within it
#' @param projection The projection genome build. Default is grch37.
#' @param min_vaf The minimum tumour VAF for a SV to be returned.
#' Default is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned.
#' Default is 40.
#' @param pass_filters If TRUE (default) only return SVs that are annotated with
#' PASS in the FILTER column. Set to FALSE to keep all variants,
#' regardless if they PASS the filters.
#' @param verbose Set to FALSE to minimize the output to console.
#' Default is TRUE. This parameter also dictates the verbose-ness of
#' any helper function internally called inside the main function.
#' @param from_cache Boolean variable for using cached results, default is TRUE.
#' If `write_to_file = TRUE`, this parameter auto-defaults to FALSE.
#' @param write_to_file Boolean statement that outputs bedpe file if TRUE,
#' default is FALSE.
#' Setting this to TRUE forces `from_cache = FALSE`.
#' @param region Specify a single region to fetch SVs anchored within
#' using the format "chrom:start-end"
#' @param these_sample_ids DEPRECATED. Use `these_samples_metadata` instead.
#' @param chromosome DEPRECATED. Use `region` instead.
#' @param qstart DEPRECATED. Use `region` instead.
#' @param qend DEPRECATED. Use `region` instead.
#' @param pairing_status DEPRECATED.
#' Subset your metadata and supply these_samples_metadata instead.

#'
#' @return A data frame in a bedpe-like format with additional
#' columns that allow filtering of high-confidence SVs.
#'
#' @import dplyr readr glue GAMBLR.helpers GAMBLR.utils
#' @export
#'
#' @examples
#' # lazily get every SV in the table with default quality filters
#' all_sv <- get_manta_sv()
#' dplyr::select(all_sv,1:14) %>% head()
#' 
#' # get all SVs for just one cohort
#' cohort_meta = suppressMessages(get_gambl_metadata()) %>% 
#'               dplyr::filter(cohort == "DLBCL_cell_lines")
#'
#' some_sv <- get_manta_sv(these_samples_metadata = cohort_meta, verbose=FALSE)
#' dplyr::select(some_sv,1:14) %>% head()
#' nrow(some_sv)
#' 
#' # get the SVs in a region around MYC
#' # WARNING: This is not the best way to find MYC SVs.
#' # Use annotate_sv on the full SV set instead.
#' myc_region_hg38 = "chr8:127710883-127761821"
#' myc_region_grch37 = "8:128723128-128774067"
#' 
#' hg38_myc_locus_sv <- get_manta_sv(region = myc_region_hg38,
#'                                 projection = "hg38",
#'                                 verbose = FALSE)
#' dplyr::select(hg38_myc_locus_sv,1:14) %>% head()
#' nrow(hg38_myc_locus_sv)
#' 
#' incorrect_myc_locus_sv <- get_manta_sv(region = myc_region_grch37,
#'                                 projection = "hg38",
#'                                 verbose = FALSE)
#' dplyr::select(incorrect_myc_locus_sv,1:14) %>% head()
#' nrow(incorrect_myc_locus_sv)
#'
#' # Despite potentially being incomplete, we can nonetheless
#' # annotate these directly for more details
#' annotated_myc_hg38 = suppressMessages(
#'          GAMBLR.utils::annotate_sv(hg38_myc_locus_sv, genome_build = "hg38")
#' )
#' head(annotated_myc_hg38)
#' table(annotated_myc_hg38$partner)
#' # The usual MYC partners are seen here
#' 
#' annotated_myc_incorrect = suppressMessages(
#'          GAMBLR.utils::annotate_sv(incorrect_myc_locus_sv, genome_build = "hg38")
#' )
#' head(annotated_myc_incorrect)
#' table(annotated_myc_incorrect$partner)
#' # The effect of specifying the wrong coordinate is evident
#' 
get_manta_sv <- function(these_samples_metadata = NULL,
                         projection = "grch37",
                         region,
                         min_vaf = 0.1,
                         min_score = 40,
                         pass_filters = TRUE,
                         verbose = TRUE,
                         from_cache = TRUE,
                         write_to_file = FALSE,
                         chromosome,
                         qstart,
                         qend,
                         these_sample_ids = NULL,
                         pairing_status) {
  if (!missing(these_sample_ids)) {
    print("parameter these_sample_ids is deprecated.
      Please use these_samples_metadata instead.")
  }
  if(missing(these_samples_metadata)){
    if(verbose){
       print("no metadata provided, fetching all samples...")
    }
    these_samples_metadata <- suppressMessages(get_gambl_metadata())
    this_meta = these_samples_metadata
  }
  if (!missing(region)) {
    region <- gsub(",", "", region)
    split_chunks <- unlist(strsplit(region, ":"))
    chromosome <- split_chunks[1]
    startend <- unlist(strsplit(split_chunks[2], "-"))
    qstart <- startend[1]
    qend <- startend[2]
  }

  if ("capture" %in% these_samples_metadata$seq_type) {
    if(verbose){
      print("dropping capture samples because manta results
      are only available for genome seq_type")
    }
    these_samples_metadata <- dplyr::filter(
      these_samples_metadata,
      seq_type == "genome"
    )
  }
  this_meta = these_samples_metadata
  if (write_to_file) {
    from_cache <- FALSE # override default automatically for nonsense combination of options
  }

  if (from_cache) {
    # get paths and check for file permissions
    output_base <- GAMBLR.helpers::check_config_and_value("project_base")
    output_file <- GAMBLR.helpers::check_config_and_value("results_merged$manta_sv$icgc_dart")
    output_file <- paste0(output_base, output_file)
    output_file <- glue::glue(output_file)

    permissions <- file.access(output_file, 4) # check read permissions

    if (permissions == -1) {
      if(verbose){
        message("No permission for unix group icgc_dart found,
        resorting to samples belonging to unix group gambl...")
      }
      output_file <- GAMBLR.helpers::check_config_and_value("results_merged$manta_sv$gambl")
      output_file <- paste0(output_base, output_file)
      output_file <- glue::glue(output_file)
    }
    if (verbose) {
      message(paste0("\nThe cached results were last updated: ", file.info(output_file)$ctime))
      message("\nReading cached results...\n")
    }

    # check for missingness of merged manta results
    if (!file.exists(output_file)) {
      print(paste("missing: ", output_file))
      message("Cannot find file locally. If working remotely, perhaps you
        forgotto load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    # read merged data
    manta_sv <- suppressMessages(read_tsv(output_file,
                                         progress = FALSE)) %>%
      dplyr::filter(
        tumour_sample_id %in% this_meta$sample_id,
        VAF_tumour >= min_vaf,
        SCORE >= min_score
      )

    if (verbose) {
      no_manta <- setdiff(this_meta$sample_id, manta_sv$tumour_sample_id)

      if (length(no_manta) > 0) {
        missing_cohorts = dplyr::filter(this_meta, sample_id %in% no_manta) %>% pull(cohort) %>% unique()
        print(paste0("No Manta SVs found for ", length(no_manta), " samples and ",length(missing_cohorts)," cohorts"))
        print(missing_cohorts)
      }
    }
  } else {
    if (write_to_file) {
      # enforce all samples in the latest metadata to be in the merge,
      # if the user decides to overwrite the cached results.
      this_meta <- suppressMessages(get_gambl_metadata()) %>%
        dplyr::filter(seq_type == "genome")
    }

    # compile the merge based on selected projection (with no filters)
    if (verbose) {
      message("\nFrom cache is set to FALSE, this function is now compiling
        a new merged results file for the selected projection...")
    }

    manta_sv <- get_manta_sv_by_samples(
      these_samples_metadata = this_meta,
      verbose = verbose,
      min_vaf = 0,
      pass_filters = FALSE,
      min_score = 0,
      projection = projection
    )

    # ensure only sample IDs in the full metadata table are kept (i.e if a
    # sample is not in the metadata table, no manta results for any such
    # sample will sneak its way into the merged results file)
    manta_sv <- manta_sv %>%
      dplyr::filter(tumour_sample_id %in% this_meta$sample_id)

    if (write_to_file) {
      # get paths and check for file permissions
      output_base <- GAMBLR.helpers::check_config_and_value("project_base")
      icgc_dart_file <- GAMBLR.helpers::check_config_and_value("results_merged$manta_sv$icgc_dart")
      icgc_dart_file <- paste0(output_base, icgc_dart_file)
      icgc_dart_file <- glue::glue(icgc_dart_file)
      icgc_dart_folder <- gsub(paste0("manta.genome--", projection, ".bedpe"), "", icgc_dart_file)

      icgc_permissions <- file.access(icgc_dart_folder, 2) # get write permission for the icgc_dart merge (all samples).

      if (icgc_permissions == 0) { # get path to gambl samples only merge, if user has acces to the icgc_dart merge.
        gambl_file <- GAMBLR.helpers::check_config_and_value("results_merged$manta_sv$gambl")
        gambl_file <- paste0(output_base, gambl_file)
        gambl_file <- glue::glue(gambl_file)

        # subset icgc_dart to only gambl samples
        gambl_samples <- this_meta %>%
          dplyr::filter(unix_group == "gambl")

        gambl_manta_sv <- manta_sv %>%
          dplyr::filter(tumour_sample_id %in% gambl_samples$sample_id)

        # write merges to file
        write_tsv(manta_sv, file = icgc_dart_file, append = FALSE)
        write_tsv(gambl_manta_sv, file = gambl_file, append = FALSE)
      } else {
        stop("You do not have sufficient permissions to write the manta merged files to disk... ")
      }
    }
  }


  # deal with chr prefixes based on the selected projection (if return is to be subset to regions...)
  if (!missing(region) || !missing(chromosome)) {
    if (projection == "grch37") {
      if (grepl("chr", chromosome)) {
        chromosome <- gsub("chr", "", chromosome)
      }
    } else if (projection == "hg38") {
      if (!grepl("chr", chromosome)) {
        chromosome <- paste0("chr", chromosome)
      }
    }

    manta_sv <- manta_sv %>%
      dplyr::filter((CHROM_A == chromosome & START_A >= qstart &
        START_A <= qend) | (CHROM_B == chromosome &
        START_B >= qstart & START_B <= qend))
  }

  if (verbose) {
    message("\nThe following VCF filters are applied;")
    message(paste0("  Minimum VAF: ", min_vaf))
    message(paste0("  Minimum Score: ", min_score))
    message(paste0("  Only keep variants passing the quality filter: ", pass_filters))
  }

  # PASS filter
  if (pass_filters) {
    manta_sv <- manta_sv %>%
      dplyr::filter(FILTER == "PASS")
  }

  # pairing status filter
  if (!missing(pairing_status)) {
    if (verbose) {
      message(paste0("  Pairing status: ", pairing_status))
    }
    manta_sv <- manta_sv %>%
      dplyr::filter(pair_status == pairing_status)
  }

  if (verbose) {
    n_variants <- nrow(manta_sv)
    unique_samples <- unique(manta_sv$tumour_sample_id)
    message(paste0(
      "\nReturning ",
      n_variants,
      " variants from ",
      length(unique_samples),
      " sample(s)"
    ))
    message("\nDone!")
  }
  #attach genome_build 
  manta_sv = dplyr::arrange(manta_sv,CHROM_A,START_A,tumour_sample_id,VAF_tumour)
  manta_sv = create_genomic_data(manta_sv, projection)
  return(manta_sv)
}
