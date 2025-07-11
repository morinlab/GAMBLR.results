#' @title Calculate Mutation Frequency By Sliding Window.
#'
#' @description Count the number of mutations in a sliding window across a region for all samples.
#'
#' @details This function is called to return the mutation frequency for a given region, either from a provided input maf data frame or from the GAMBL maf data.
#' @details Regions are specified with the `region`parameter.
#' @details Alternatively, the region of interest can also be specified by calling the function with `chromosome`, `start_pos`, and `end_pos` parameters.
#' @details This function operates on a single region. To return a matrix of sliding window counts over multiple regions, see [GAMBLR.results::calc_mutation_frequency_bin_region]Use .
#'
#' @param region A string describing a genomic region in the "chrom:start-end" format. The region must be specifed in this format OR as separate chromosome, start_pos, end_pos arguments.
#' @param chromosome Chromosome name in region.
#' @param start_pos Start coordinate of region.
#' @param end_pos End coordinate of region.
#' @param these_samples_metadata Optional data frame containing a sample_id column. If not providing a maf file, seq_type is also a required column.
#' @param these_sample_ids Optional vector of sample IDs. Output will be subset to IDs present in this vector.
#' @param maf_data Optional maf data frame. Will be subset to rows where Tumor_Sample_Barcode matches provided sample IDs or metadata table. If not provided, maf data will be obtained with get_ssm_by_regions().
#' @param projection Specify which genome build to use. Required.
#' @param slide_by Slide size for sliding window. Default 100.
#' @param window_size Size of sliding window. Default 1000.
#' @param return_format Return format of mutations. Accepted inputs are "long" and "wide". Long returns a data frame of one sample ID/window per row. Wide returns a matrix with one sample ID per row and one window per column. Using the "wide" format will retain all samples and windows regardless of the drop_unmutated or min_count_per_bin parameters.
#' @param min_count_per_bin Minimum counts per bin, default is 0. Setting this greater than 0 will drop unmutated windows only when return_format is long.
#' @param return_count Boolean statement to return mutation count per window (TRUE) or binary mutated/unmutated status (FALSE). Default is TRUE.
#' @param drop_unmutated Boolean for whether to drop windows with 0 mutations. Only effective with "long" return format.
#'
#' @return Either a matrix or a long tidy table of counts per window.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' meta = suppressMessages(get_gambl_metadata()) %>% 
#'                         dplyr::filter(pathology=="MCL")
#' 
#' mut_freq = calc_mutation_frequency_bin_region(these_samples_metadata = meta,
#'                                               region = "11:69455000-69459900",
#'                                               slide_by = 10,
#'                                               window_size = 10000)
#' head(mut_freq)
#' 
#' \dontrun{
#' # This will fail because the chromosome naming doesn't match the default projection 
#' misguided_attempt = calc_mutation_frequency_bin_region(these_samples_metadata = meta,
#'                                                          region = "chr11:69455000-69459900",
#'                                                          slide_by = 10,
#'                                                          window_size = 10000) 
#' # This will work!
#' mut_freq = calc_mutation_frequency_bin_region(these_samples_metadata = meta,
#'                                                          region = "chr11:69455000-69459900",
#'                                                          slide_by = 10,
#'                                                          window_size = 10000,projection="hg38")
#' head(mut_freq)
#' }
calc_mutation_frequency_bin_region <- function(region,
                                          chromosome,
                                          start_pos,
                                          end_pos,
                                          these_samples_metadata,
                                          these_sample_ids = NULL,
                                          maf_data = NULL,
                                          projection = "grch37",
                                          slide_by = 100,
                                          window_size = 1000,
                                          return_format = "long",
                                          min_count_per_bin = 0,
                                          return_count = TRUE,
                                          drop_unmutated = FALSE
                                          ) {
  if(missing(these_samples_metadata)){
    stop("metadata must be provided via the these_samples_metadata argument")
  }
  these_samples_metadata = dplyr::filter(these_samples_metadata, seq_type %in% c("genome", "capture"))
  if(dim(these_samples_metadata)[1] == 0){
    stop("these_samples_metadata must have samples with seq_type genome or capture")
  }
  # Create objects to describe region both as string and individual objects
  try(if (missing(region) & missing(chromosome)) {
    stop("No region information provided. Please provide a region as a string in the chrom:start-end format, or as individual arguments. ")
  })

  if ((drop_unmutated | min_count_per_bin > 0) & return_format == "wide") {
    message("To return a wide table, all samples and windows must be kept. Ignoring drop_unmutated and min_count_per_bin arguments. ")
  }

  if (missing(region)) {
    region <- paste0(
      chromosome, ":", start_pos, "-",
      end_pos
    )
  } else {
    chunks <- region_to_chunks(region)
    chromosome <- chunks$chromosome
    start_pos <- as.numeric(chunks$start)
    end_pos <- as.numeric(chunks$end)
  }
  
  metadata <- these_samples_metadata
  if(!missing(these_sample_ids)){
    stop("deprecated agrument these_sample_ids was provided. Please use these_samples_metadata instead")
  }
  # these_sample_ids <- pull(metadata,sample_id) %>% unique


  if (
    (str_detect(chromosome, "chr") & projection == "grch37") |
      (!str_detect(chromosome, "chr") & projection == "hg38")
  ) {
    stop("chr prefixing status of region and specified projection don't match. ")
  }


  # Check region size and compare to max region size
  # Is this really needed?
  max_region <- 5e+06

  region_size <- end_pos - start_pos
  if (region_size < max_region) {
    message(paste(
      "processing bins of size", window_size,
      "across", region_size, "bp region"
    ))
  } else {
    message(paste("CAUTION!\n", region_size, "exceeds maximum size recommended by this function."))
  }

  # Split region into windows
  windows <- data.frame(
    chrom = chromosome,
    window_start = seq(start_pos, end_pos, by = slide_by)
  ) %>%
    dplyr::mutate(window_end = window_start + window_size - 1) %>%
    dplyr::select(chrom, window_start, window_end)

  # Option to return full region count instead of sliding window
  if (window_size == 0) {
    windows <- data.frame(
      chrom = chromosome,
      window_start = start_pos,
      window_end = end_pos
    )
  }

  # Obtain SSM coordinates from GAMBL if no maf_data was provided
  if (is.null(maf_data)) {
    try(
      if (!"seq_type" %in% colnames(metadata)) {
        stop("seq_type must be present in metdata for compatibility with get_ssm_by_sample")
      }
    )
    message("Using GAMBLR.results::get_ssm_by_region...")
    seq_type_sample_ids = list() # items are dataframes
    for(a_seq_type in unique(these_samples_metadata$seq_type)){
        seq_type_sample_ids[[a_seq_type]]=dplyr::filter(these_samples_metadata,seq_type==a_seq_type)
    }
    seq_type_region_ssm <- list()
    for (a_seq_type in names(seq_type_sample_ids)) {
      seq_type_region_ssm[[a_seq_type]] <- GAMBLR.results:::get_ssm_by_region(
        region = region,
        projection = projection,
        streamlined = FALSE,
        these_samples_metadata = seq_type_sample_ids[[a_seq_type]]
      ) %>%
        dplyr::mutate(end = Start_Position + 1) %>%
        dplyr::select(
          chrom = Chromosome,
          start = Start_Position,
          end,
          sample_id = Tumor_Sample_Barcode
        )
    }

    region_ssm <-  do.call("rbind", seq_type_region_ssm)

  } else {
    #  Subset provided maf to specified region
    message("Using provided maf...")
    maf.dt <- data.table(maf_data)
    region_bed <- data.table(
      "Chromosome" = as.character(chromosome),
      "Start_Position" = as.numeric(start_pos),
      "End_Position" = as.numeric(end_pos)
    )
    setkey(region_bed)
    region_ssm <- foverlaps(maf.dt, region_bed) %>%
      dplyr::filter(!is.na(Start_Position)) %>%
      dplyr::mutate(end = i.Start_Position - 1) %>%
      dplyr::select(
        chrom = Chromosome,
        start = i.Start_Position,
        end,
        sample_id = Tumor_Sample_Barcode
      ) %>%
      dplyr::mutate(mutated = 1)

    region_ssm <- data.frame(metadata) %>%
      dplyr::select(sample_id) %>%
      dplyr::left_join(region_ssm) %>%
      dplyr::filter(!is.na(mutated))
  }

  # Check if the region is empty.
  # If yes return NULL so that running this function with lapply will allow bind_rows to run on the output.
  if (nrow(region_ssm) == 0 & (drop_unmutated | min_count_per_bin > 0)) {
    message(paste0("No mutations found in region ", region, " for this sample set. "))
    return(NULL)
  }

  # Count mutations per window
  windows_tallied <- dplyr::inner_join(
    windows,
    region_ssm,
    by = "chrom"
  ) %>%
    dplyr::filter(
      start >= window_start,
      start <= window_end
    ) %>%
    dplyr::group_by(
      sample_id,
      window_start
    ) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::full_join(select(metadata, sample_id)) %>%
    dplyr::arrange(sample_id) %>%
    dplyr::full_join(select(windows, window_start)) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = window_start,
      values_from = n,
      values_fill = 0
    ) %>%
    dplyr::select(-matches("^NA$")) %>%
    tidyr::pivot_longer(
      -c(sample_id),
      names_to = "window_start",
      values_to = "n"
    ) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(sample_id))

  # Remove unmutated windows if requested
  if (drop_unmutated | min_count_per_bin > 0) {
    windows_tallied <- windows_tallied %>%
      dplyr::filter(n >= min_count_per_bin)
    if (drop_unmutated & min_count_per_bin == 0) {
      windows_tallied %>%
        dplyr::filter(n > 0)
    }
  }

  # Create requested data output format
  if (return_count) {
    # Return table of mutation counts per bin
    windows_tallied_final <- mutate(
      windows_tallied,
      bin = paste0(chromosome, "_", window_start)
    ) %>%
      dplyr::mutate(mutation_count = n) %>%
      dplyr::select(
        sample_id,
        bin,
        mutation_count
      )
  } else {
    # Return table of binary mutated/unmutated status per bin
    windows_tallied_final <- mutate(
      windows_tallied,
      bin = paste0(chromosome, "_", window_start)
    ) %>%
      dplyr::mutate(mutated = ifelse(n > 0, 1, 0)) %>%
      dplyr::select(
        sample_id,
        bin,
        mutated
      )
  }

  if (return_format == "wide") {
    widened <- windows_tallied_final %>%
      tidyr::pivot_wider(
        names_from = bin,
        values_from = matches("mutat"),
        values_fill = 0
      )
    return(widened)
  } else {
    return(windows_tallied_final)
  }
}
