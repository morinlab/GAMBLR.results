

#' @title Get Coding SSM.
#'
#' @description Retrieve all coding SSMs from one seq_type in GAMBL
#' in MAF-like format.
#'
#' @details Effectively retrieve simple somatic mutations (SSM) results for
#' either capture or genome seq_type (but not both at once).
#' The resulting data frame will be a maf_data object, which tracks the
#' genome build (projection) for the variants and will have a
#' maf_seq_type column that tracks the origin seq_type each variant.
#' In most cases, users should be using the related function that
#' is able to obtain SSMs across both genome and capture seq_type:
#' [GAMBLR.results::get_all_coding_ssm] 
#' Is this function not what you are looking for? Try one of:
#' [GAMBLR.results::get_coding_ssm_status],
#' [GAMBLR.results::get_ssm_by_patients],
#' [GAMBLR.results::get_ssm_by_sample],
#' [GAMBLR.results::get_ssm_by_samples],
#' [GAMBLR.results::get_ssm_by_region],
#' [GAMBLR.results::get_ssm_by_regions]
#'
#' @param these_samples_metadata Optional (but highly recommended) metadata
#' table to tell the function how to subset the data. Only sample_id in
#' this table with the matching seq_type will be in the output. Not all
#' samples may be in the output, though, if they are missing SSM results
#' or had no mutations detected.
#' @param force_unmatched_samples Optional argument for forcing unmatched
#' samples, using [GAMBLR.results::get_ssm_by_samples].
#' @param projection Reference genome build for the coordinates in the MAF file.
#' The default is grch37.
#' @param this_seq_type The seq_type you want SSMs from, default is genome.
#' @param basic_columns Basic columns refers to the first 45 standard MAF
#' columns. Set this to FALSE if you want all the available columns instead.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify
#' which columns to be returned within the MAF. This parameter can either
#' be a vector of indexes (integer) or a vector of characters (matching
#' columns in MAF).
#' @param augmented Set to FALSE if you instead want the original MAF from
#' each sample for multi-sample patients instead of the augmented MAF
#' @param min_read_support Only returns variants with at least this many
#' reads in t_alt_count (for cleaning up augmented MAFs)
#' @param include_silent Logical indicating whether to include
#' silent mutations in coding regions (i.e. synonymous). Default is TRUE.
#' @param verbose Controls the "verboseness" of this function (and
#' internally called helpers).
#' @param from_flatfile Deprecated. Ignored
#' @param engine Deprecated. Ignored
#' @param groups Deprecated. Use these_samples_metadata instead.
#' @param these_sample_ids Deprecated. Use these_samples_metadata instead.
#' @param limit_cohort Deprecated. Use these_samples_metadata instead.
#' @param exclude_cohort  Deprecated. Use these_samples_metadata instead.
#' @param limit_pathology Deprecated. Use these_samples_metadata instead.
#' @param limit_samples Deprecated. Use these_samples_metadata instead.
#'
#' @return A data frame containing all the MAF data columns
#' (one row per mutation).
#'
#' @import dplyr tidyr RMariaDB DBI glue GAMBLR.helpers GAMBLR.utils
#' @export
#'
#' @examples
#' 
#'   #basic usage (defaults to genome seq_type)
#'   maf_genome = get_coding_ssm()
#' 
#'   nrow(maf_genome)
#' 
#'   dplyr::select(maf_genome,1,4,5,6,9,maf_seq_type)
#' 
#'   maf_exome_hg38 = get_coding_ssm(this_seq_type = "capture",
#'                                   projection="hg38") 
#' 
#'   dplyr::select(maf_exome_hg38,1,4,5,6,9,maf_seq_type)
#' 
get_coding_ssm = function(these_samples_metadata = NULL,
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
                          verbose = FALSE,
                          limit_cohort,
                          exclude_cohort,
                          limit_pathology,
                          limit_samples,
                          from_flatfile,
                          these_sample_ids){
  

  if(any(!missing(groups),!missing(limit_cohort),!missing(exclude_cohort),!missing(limit_pathology),!missing(limit_samples),!missing(these_sample_ids))){
    stop("limit_samples, limit_cohort, exclude_cohort and limit_pathology, these_sample_ids are deprecated. Use `these_samples_metadata` instead")
  }
  
  remote_session = check_remote_configuration()
  
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }

  

  #get file path for augmented maf
  if(augmented){
    seq_type = this_seq_type #glue needs to have this variable
    maf_template = check_config_and_value("results_flatfiles$ssm$template$cds$augmented")
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(check_config_and_value("project_base"), maf_path)
  }else{
    seq_type = this_seq_type #glue needs to have this variable
    maf_template = check_config_and_value("results_flatfiles$ssm$template$cds$deblacklisted")
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(check_config_and_value("project_base"), maf_path)
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
    muts = suppressMessages(fread_maf(full_maf_path,select_cols=select_cols)) %>%
        dplyr::filter(Variant_Classification %in% coding_class) %>%
        as.data.frame()
  }else{
    muts = suppressMessages(fread_maf(full_maf_path)) %>%
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


  if(!missing(these_samples_metadata)){
      these_samples_metadata = dplyr::filter(these_samples_metadata,seq_type == this_seq_type)
      muts = dplyr::filter(muts,Tumor_Sample_Barcode %in% these_samples_metadata$sample_id)
      mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
      if(verbose){
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
  
  muts = mutate(muts,maf_seq_type = this_seq_type)
  muts = GAMBLR.utils::create_maf_data(muts,projection)
  return(muts)
}
