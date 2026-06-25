#' @title Maf To Custom Track.
#'
#' @description Convert mutations into a UCSC custom track file
#'
#' @details This function takes a set of mutations as maf_data and converts it
#' to a UCSC Genome Browser ready BED (or bigbed/biglolly) file complete with
#' the required header. Upload the resulting file to
#' [UCSC genome browser](<https://genome.ucsc.edu/cgi-bin/hgCustom>)
#' to view your data as a custom track.
#' Optional parameters available for further customization of the returned file.
#' For more information, refer to the parameter descriptions
#' and function examples.
#'
#' @param maf_data MAF data obtained from of the `get_ssm` family of functions.
#' @param these_samples_metadata A metadata table to subset the samples of interest from
#'   the input `maf_data`. If NULL (the default), all samples in `maf_data` are kept.
#' @param this_seq_type The seq type used to filter the samples if `these_samples_metadata` 
#'   is not provided. Should match what was used to get `maf_data`. Default is "genome".
#' @param output_file Name for your new bed file that can be uploaded as a custom track to UCSC.
#' @param as_bigbed Boolean parameter controlling the format of the created track file.
#'   Default is FALSE.
#' @param as_biglolly Boolean parameter controlling the format of the created track file.
#'   Default is FALSE (i.e a BED file will be returned).
#' @param colour_column Set the colouring the SSMs in the track. Possible options are lymphgen
#'   (default), pathology, genome_build (as in `these_samples_metadata`), and mutation
#'   (corresponds to MAF Variant_Classification).
#' @param track_name Track name to use in the header if output is not bigBed or bigLolly. 
#'    Default is "GAMBL mutations"
#' @param track_description Track description to use in the header if output is not 
#'    bigBed or bigLolly. Default is "GAMBL mutations"
#' @param verbose Default is FALSE.
#' @param padding_size Optional parameter specifying the padding size in the
#'   returned file, default is 0.
#' @param projection Specify which genome build to use. Possible values are "grch37" (default)
#'    or "hg38". This parameter has an effect only when `as_bigbed` or `as_biglolly` is TRUE.
#' @param bedToBigBed_path Path to your local `bedToBigBed` UCSC tool or the string
#'   `"config"` (default). If set to `"config"`, `GAMBLR.helpers::check_config_value`
#'   is called internally and the `bedToBigBed` path is obtained from the `config.yml`
#'   file saved in the current working directory. This parameter is ignored if both
#'   `as_bigbed` and `as_biglolly` is set to `FALSE`.
#' @param these_sample_ids DEPRECATED
#' @return Nothing.
#'
#' @import tidyr dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' # using grch37 coordinates
#' myc_grch37 <- GAMBLR.utils::create_bed_data(
#'                 GAMBLR.data::grch37_lymphoma_genes_bed
#'               ) %>%
#'               dplyr::filter(name == "MYC")
#'
#' print(myc_grch37)
#' genomes <- get_gambl_metadata() %>% dplyr::filter(seq_type == "genome")
#' genome_maf <- get_ssm_by_regions(regions_bed = myc_grch37,
#'                              these_samples_metadata = genomes,
#'                              streamlined = FALSE,
#'                              basic_columns = TRUE)
#'
#' # myc_hg19.bed will be created in your working directory
#' maf_to_custom_track(maf_data = genome_maf, output_file = "myc_genome_hg19.bed")
#'
#' #lazy/concise way:
#' my_region = "8:128747680-128753674"
#'
#' captures <- get_gambl_metadata() %>% dplyr::filter(seq_type == "capture")
#' capture_maf <- get_ssm_by_regions(regions_list = my_region,
#'                              these_samples_metadata = captures,
#'                              projection = "grch37",
#'                              streamlined = FALSE,
#'                              basic_columns = TRUE)
#' maf_to_custom_track(maf_data = capture_maf, this_seq_type = "capture", output_file = "myc_capture_hg19.bed")
#'
maf_to_custom_track <- function(maf_data,
                                these_samples_metadata = NULL,
                                this_seq_type = "genome",
                                output_file,
                                as_bigbed = FALSE,
                                as_biglolly = FALSE,
                                colour_column = "lymphgen",
                                track_name = "GAMBL mutations",
                                track_description = "GAMBL mutations",
                                verbose = FALSE,
                                padding_size = 0,
                                projection = "grch37",
                                bedToBigBed_path = "config",
                                these_sample_ids = NULL) {
  # Handle deprecated parameters
  if (!is.null(these_sample_ids)) {
    stop("Parameter `these_sample_ids` is deprecated and will be ignored.
    Please use `these_samples_metadata` instead.")
  }
  if (missing(these_samples_metadata)) {
    these_samples_metadata <- get_gambl_metadata() %>%
      dplyr::filter(seq_type == this_seq_type)
  }
  # check some provided parameter
  if (as_bigbed) {
    if (bedToBigBed_path == "config") {
      bedToBigBed_path <- tryCatch(
        check_config_and_value(config::get("dependencies")$bedToBigBed),
        error = function(e) {
          k <- paste0("You set bedToBigBed_path parameter to \"config\". However...\n", e)
          stop(k, call. = FALSE)
        }
      )
    } else {
      stopifnot(
        "`bedToBigBed_path` points to a non-existent file." =
          file.exists(bedToBigBed_path)
      )
    }
  }

  # subset input maf according to metadata samples
  maf_data <- dplyr::filter(maf_data, Tumor_Sample_Barcode %in% these_samples_metadata$sample_id)

  # reduce to a bed-like format
  if(colour_column %in% "mutation"){
    maf_data <- dplyr::select(maf_data, Chromosome, Start_Position, End_Position, Variant_Classification)
    colnames(maf_data) <- c("chrom", "start", "end", "group")
  }else{
    maf_data <- dplyr::select(maf_data, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)
    colnames(maf_data) <- c("chrom", "start", "end", "sample_id")
  }

  maf_data <- dplyr::mutate(maf_data, end = end + padding_size)
  if (!any(grepl("chr", maf_data[, 1]))) {
    # add chr prefix
    maf_data[, 1] <- unlist(lapply(maf_data[, 1], function(x) {
      paste0("chr", x)
    }))
  }

  stopifnot("`colour_column` must be one of \"lymphgen\", \"pathology\", \"genome_build\", or \"mutation\"" =
    colour_column %in% c("lymphgen", "pathology", "genome_build", "mutation"))

  colour_cols <- GAMBLR.helpers::get_gambl_colours(colour_column, verbose = verbose)

  colour_df <- data.frame(group = names(colour_cols), colour = colour_cols)
  rgb_df <- data.frame(t(col2rgb(colour_cols))) %>%
    dplyr::mutate(group = names(colour_cols), hex = unname(colour_cols)) %>%
    tidyr::unite(col = "rgb", red, green, blue, sep = ",")
  if (verbose) {
    print(rgb_df)
  }

  # colours associated with samples if colour_column is not mutation
  # otherwise they are associated with the variants
  if(!colour_column %in% "mutation"){
    meta <- dplyr::select(these_samples_metadata, sample_id, all_of(colour_column))
    colnames(meta)[2] <- "group"
    samples_coloured <- left_join(meta, rgb_df) %>%
      dplyr::mutate(group = ifelse(is.na(group), "NotAvailable", group),
            rgb = ifelse(is.na(rgb), "0,0,0", rgb))
    if (verbose) {
      print(samples_coloured)
    }
    maf_bed <- maf_data %>%
      dplyr::mutate(score = 0, strand = ".", thickStart = start - 1, start = thickStart, thickEnd = end)
    if (verbose) {
      print(head(maf_bed))
    }
    maf_coloured <- left_join(maf_bed, samples_coloured, by = "sample_id") %>%
      dplyr::select(chrom, start, end, group, score, strand, thickStart, thickEnd, rgb)
  }else{
    maf_coloured <- left_join(maf_data, rgb_df) %>%
      dplyr::mutate(group = ifelse(is.na(group), "Unknown", group),
            rgb = ifelse(is.na(rgb), "0,0,0", rgb)) %>%
      dplyr::mutate(score = 0, strand = ".", thickStart = start - 1, start = thickStart, thickEnd = end) %>%
      dplyr::select(chrom, start, end, group, score, strand, thickStart, thickEnd, rgb)
  }

  maf_summary <- group_by(maf_coloured, rgb) %>% tally()
  if (verbose) {
    print(maf_summary)
    print(head(maf_coloured))
  }

  if (as_bigbed | as_biglolly) {
    if (grepl(pattern = ".bb$", x = output_file)) {
      # temp file will be .bed
      temp_bed <- tempfile(pattern = "bed_")
    } else {
      stop("please provide an output file name ending in .bb
      to create a bigBed file")
    }

    maf_coloured <- maf_coloured %>%
      arrange(chrom, start)

    # create temp file chrom.sizes
    if (projection == "grch37") {
      chr_arms <- GAMBLR.data::chromosome_arms_grch37 %>%
        dplyr::mutate(chromosome = paste0("chr", chromosome))
    } else if (projection == "hg38") {
      chr_arms <- GAMBLR.data::chromosome_arms_hg38
    } else {
      stop("projection parameter must be \"grch37\" or \"hg38\".")
    }
    chr_sizes <- chr_arms %>%
      dplyr::filter(arm == "q") %>%
      dplyr::select(chromosome, end) %>%
      dplyr::rename(size = end)
    temp_chr_sizes <- tempfile(pattern = "chrom.sizes_")
    write_tsv(chr_sizes, file = temp_chr_sizes, col_names = FALSE)

    if (as_biglolly) {
      # currently the same code is run either way but this may change so I've separated this until we settle on format
      # TO DO: collapse based on hot spot definition and update column 4 (score) based on recurrence
      # needs to have size column
      maf_score_options <- factor(maf_coloured$rgb)
      maf_coloured$score <- as.numeric(maf_score_options)

      # determine frequency of each event per group to assign the size
      maf_coloured <- group_by(maf_coloured, start, rgb) %>% mutate(size = n())

      write_tsv(maf_coloured, file = temp_bed, col_names = FALSE)
      # conversion:
      autosql_file <- "/Users/rmorin/git/LLMPP/resources/reference/ucsc/bigLollyExample3.as"

      bigbed_conversion <- paste(
        bedToBigBed_path, "-as=", autosql_file, "-tab -type=bed9+1", temp_bed,
        temp_chr_sizes, output_file
      )
      print(bigbed_conversion)
      system(bigbed_conversion)
    } else {
      write_tsv(maf_coloured, file = temp_bed, col_names = F)
      # conversion:
      bigbed_conversion <- paste(bedToBigBed_path, "-tab -type=bed9", temp_bed, temp_chr_sizes, output_file)

      system(bigbed_conversion)
    }
    unlink(c(temp_bed, temp_chr_sizes))
  } else {
    header_ucsc <- paste0('track name="', track_name, '" description="', track_description, '" visibility=2 itemRgb="On"\n')
    cat(header_ucsc, file = output_file)
    write_tsv(maf_coloured, file = output_file, col_names = F, append = TRUE)
  }
}
