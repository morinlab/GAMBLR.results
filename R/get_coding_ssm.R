
#' Create MAF Data
#'
#' This function creates MAF (Mutation Annotation Format) data from the given input.
#'
#' @param maf_df A data frame containing the MAF data.
#' @param genome_build A string specifying the genome build ("grch37" or "hg38").
#' @return A data frame with class attributes for MAF data.
#' @export
create_maf_data <- function(maf_df, genome_build) {
  if (!inherits(maf_df, "data.frame")) stop("data must be a data frame")
  if (!genome_build %in% c("grch37", "hg38")) stop("Invalid genome build")
  
  structure(maf_df,
            class = c("maf_data", "genomic_data", class(maf_df)),  #  "genomic_data" for generic methods
            genome_build = genome_build)
}

#' Get Genome Build
#'
#' This function retrieves the genome build attribute from the data.
#'
#' @param data A data frame with genome build attribute.
#' @return A string specifying the genome build.
#' @export
get_genome_build <- function(data) {
  attr(data, "genome_build")
}

#' Preserve Genomic Attributes
#'
#' This function preserves the genomic attributes and class after dplyr operations.
#'
#' @param new_data A data frame resulting from dplyr operations.
#' @param old_data The original data frame with genomic attributes.
#' @return A data frame with preserved genomic attributes.
#' @export
preserve_genomic_attributes <- function(new_data, old_data) {
  # Preserve the genome_build attribute
  attr(new_data, "genome_build") <- attr(old_data, "genome_build")
  
  # Combine the new data’s classes with the genomic classes
  new_data_classes <- class(new_data)
  # Ensure the genomic classes are at the front
  new_classes <- unique(c("maf_data", "genomic_data", new_data_classes))
  class(new_data) <- new_classes
  
  return(new_data)
}

# S3 methods for genomic_data class
#' @export
mutate.genomic_data <- function(.data, ...) {
  new_data <- dplyr::mutate(as.data.frame(.data), ...)  
  preserve_genomic_attributes(new_data, .data)
}
#' @export
filter.genomic_data <- function(.data, ...) {
  new_data <- dplyr::filter(as.data.frame(.data), ...)  
  preserve_genomic_attributes(new_data, .data)
}
#' @export
select.genomic_data <- function(.data, ...) {
  new_data <- dplyr::select(as.data.frame(.data), ...)  
  preserve_genomic_attributes(new_data, .data)
}
#' @export
rename.genomic_data <- function(.data, ...) {
  new_data <- dplyr::rename(as.data.frame(.data), ...)  
  preserve_genomic_attributes(new_data, .data)
}
#' @export
arrange.genomic_data <- function(.data, ...) {
  new_data <- dplyr::arrange(as.data.frame(.data), ...)  
  preserve_genomic_attributes(new_data, .data)
}
#' @export
group_by.genomic_data <- function(.data, ..., .add = FALSE) {
  new_data <- dplyr::group_by(as.data.frame(.data), ..., .add = .add)  
  preserve_genomic_attributes(new_data, .data)
}
#' @export
ungroup.genomic_data <- function(x, ...) {
  new_data <- dplyr::ungroup(as.data.frame(x), ...)
  preserve_genomic_attributes(new_data, x)
}
# Merger function to preserve metadata and protect against accidental mixing across genome_builds

#' Bind maf or other genomic data together
#'
#' @description Combine multiple maf_data objects and retain metadata such as genome_build.
#' This function will not allow you to combine maf_data objects that have different genome_build values.
#' An error will also be thrown if the same sample id is found in more than one of the inputs (if check_id is TRUE).
#'
#' @param ... All maf_data or seg_data objects to be combined.
#' @param check_id Logical. If TRUE (the default), the function will check for the presence of the expected ID column
#'        and for duplicate sample IDs across the inputs. Set to FALSE to skip this check.
#'
#' @return data.frame with combined data and preserved genome_build metadata.
#' @export
#'
#' @examples
#'
#' merged_maf = bind_genomic_data(maf1, maf2,check_id=FALSE)
#'
bind_genomic_data <- function(..., check_id = TRUE) {
  
  in_list <- list(...)
  
  if ("maf_data" %in% class(in_list[[1]])) {
    # MAF format, ID column is Tumor_Sample_Barcode
    id_col <- "Tumor_Sample_Barcode"
  } else if ("seg_data" %in% class(in_list[[1]])) {
    # SEG format, ID column is ID
    id_col <- "ID"
  } else {
    stop(paste("Unsure how to merge:", class(in_list[[1]])))
  }
  
  # Ensure all inputs are either maf_data or seg_data objects
  if (!all(sapply(in_list, inherits, "maf_data")) && !all(sapply(in_list, inherits, "seg_data"))) {
    stop("All inputs must be maf_data objects or seg_data objects.")
  }
  
  # Extract genome builds
  genome_builds <- unique(sapply(in_list, get_genome_build))
  
  if (length(genome_builds) > 1) {
    stop("Cannot bind seg_data or maf_data objects with different genome builds: ", 
         paste(genome_builds, collapse = ", "))
  }
  
  # If check_id is TRUE, verify that the expected ID column exists and that IDs are unique.
  if (check_id) {
    # Collect unique sample IDs from each dataset
    id_sets <- lapply(in_list, function(df) {
      if (!(id_col %in% colnames(df))) {
        stop("ID column '", id_col, "' not found in input data.")
      }
      unique(df[[id_col]])
    })
    
    # Flatten the list and count occurrences of each ID
    all_ids <- unlist(id_sets)
    duplicate_ids <- names(table(all_ids)[table(all_ids) > 1])
    
    # If any ID is found in multiple datasets, throw an error
    if (length(duplicate_ids) > 0) {
      stop("Duplicate IDs found in multiple input data frames: ", paste(duplicate_ids, collapse = ", "))
    }
  }
  
  combined <- dplyr::bind_rows(in_list)
  attr(combined, "genome_build") <- genome_builds[1]  # Assign the common genome build
  
  if (!"maf_data" %in% class(combined)) {
    class(combined) <- c("maf_data", "genomic_data", class(combined))  # Preserve class
  }
  
  return(combined)
}


#' @title Get Coding SSM.
#'
#' @description Retrieve all coding SSMs from one seq_type in GAMBL in MAF-like format.
#'
#' @details Effectively retrieve coding SSM calls. Multiple filtering parameters are available for this function.
#' For more information on how to implement the filtering parameters, refer to the parameter descriptions as well as examples in the vignettes.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_coding_ssm_status],
#' [GAMBLR.results::get_ssm_by_patients], [GAMBLR.results::get_ssm_by_sample], [GAMBLR.results::get_ssm_by_samples], [GAMBLR.results::get_ssm_by_region], [GAMBLR.results::get_ssm_by_regions]
#'

#' @param these_samples_metadata Supply a metadata table to auto-subset the data to samples in that table before returning.
#' @param these_sample_ids Optional, restrict the returned maf to a set of sample IDs (this parameter is replacing the `limit_samples` parameter).
#' @param force_unmatched_samples Optional argument for forcing unmatched samples, using [GAMBLR.results::get_ssm_by_samples].
#' @param projection Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#' @param this_seq_type The seq_type you want back, default is genome.
#' @param basic_columns Set to FALSE to override the default behavior of returning only the first 45 columns of MAF data.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters (matching columns in MAF).
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs)
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is TRUE.
#' @param verbose Controls the "verboseness" of this function (and internally called helpers).
#' @param from_flatfile Deprecated. Ignored
#' @param engine Deprecated. Ignored
#' @param groups Deprecated. Use either these_samples_metadata or these_sample_ids instead 
#' @param limit_cohort Deprecated. Use either these_samples_metadata or these_sample_ids instead 
#' @param exclude_cohort  Deprecated. Use either these_samples_metadata or these_sample_ids instead 
#' @param limit_pathology Deprecated. Use either these_samples_metadata or these_sample_ids instead 
#' @param limit_samples Deprecated. Use either these_samples_metadata or these_sample_ids instead 
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#'
#' @import dplyr tidyr RMariaDB DBI glue GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#'   #basic usage
#'   maf_data = get_coding_ssm(limit_cohort = c("BL_ICGC"))
#'
#'   maf_data = get_coding_ssm(limit_samples = "HTMCP-01-06-00485-01A-01D")
#' }
get_coding_ssm = function(
                          these_sample_ids = NULL,
                          these_samples_metadata = NULL,
                          force_unmatched_samples,
                          projection = "grch37",
                          this_seq_type = "genome",
                          basic_columns = TRUE,
                          maf_cols = NULL,
                          augmented = TRUE,
                          min_read_support = 3,
                          groups = c("gambl", "icgc_dart"),
                          include_silent = TRUE,
                          engine,
                          verbose = TRUE,
                          limit_cohort,
                          exclude_cohort,
                          limit_pathology,
                          limit_samples,
                          from_flatfile){
  
  if(!missing(limit_samples)){
    
  }
  if(any(!missing(groups),!missing(limit_cohort),!missing(exclude_cohort),!missing(limit_pathology),!missing(limit_samples))){
    stop("limit_samples, limit_cohort, exclude_cohort and limit_pathology are deprecated, use `these_sample_ids` instead, or use a metadata table already subset to the samples of interest with `these_samples_metadata`")
  }

  remote_session = check_remote_configuration()
  
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }

  

  #get file path for augmented maf
  if(augmented){
    seq_type = this_seq_type #glue needs to have this variable
    maf_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$ssm$template$cds$augmented)
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), maf_path)
  }else{
    seq_type = this_seq_type #glue needs to have this variable
    maf_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$ssm$template$cds$deblacklisted)
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), maf_path)
  }

  #read file

  if(verbose){
    message(paste("reading from:", full_maf_path)) 
  }

  #check for missingness
  if(!file.exists(full_maf_path)){
    print(paste("missing: ", full_maf_path))
    message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
    message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
  }


  if(basic_columns){
    #subset to basic columns during read to save time and memory with lazy loading (in theory)
    select_cols = c(1:45)
    muts = fread_maf(full_maf_path,select_cols=select_cols) %>%
        dplyr::filter(Variant_Classification %in% coding_class) %>%
        as.data.frame()
  }else{
    muts = fread_maf(full_maf_path) %>%
        dplyr::filter(Variant_Classification %in% coding_class) %>%
        as.data.frame()
  }


  if(verbose){
    mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
    message(paste("mutations from", mutated_samples, "samples")) 
  }

  #if augmented maf selected, drop variants with low read support (default is 3)
  if(augmented){
    muts = dplyr::filter(muts, t_alt_count >= min_read_support)
  }

  #filter maf on selected sample ids
  if(!missing(these_sample_ids)){
    muts = muts %>%
      dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids)
  }

  if(verbose){
    if(!missing(these_samples_metadata)){
      muts = dplyr::filter(muts,Tumor_Sample_Barcode %in% these_samples_metadata$sample_id)
      mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
      message(paste("after linking with metadata, we have mutations from exactly", mutated_samples, "samples")) 
    }
    
  }

  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    muts = dplyr::select(muts, all_of(maf_cols))
  }

  #drop rows for these samples so we can swap in the force_unmatched outputs instead
  if(!missing(force_unmatched_samples)){
    muts = muts %>%
      dplyr::filter(!Tumor_Sample_Barcode %in% force_unmatched_samples)

    if(verbose){
      nsamp = length(force_unmatched_samples)
      message(paste("dropping variants from", nsamp, "samples and replacing with force_unmatched outputs")) 
    }

    #get replacements using get_ssm_by_samples
    fu_muts = get_ssm_by_samples(these_sample_ids = force_unmatched_samples)
    muts = bind_rows(muts, fu_muts)
  }
  
  muts = create_maf_data(muts,projection)
  # use S3-safe version of dplyr function
  muts = mutate.genomic_data(muts,maf_seq_type = this_seq_type)
  return(muts)
}
