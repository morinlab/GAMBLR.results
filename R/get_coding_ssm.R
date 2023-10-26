#' @title Get Coding SSM.
#'
#' @description Retrieve all coding SSMs from the GAMBL database in MAF-like format.
#'
#' @details Effectively retrieve coding SSM calls. Multiple filtering parameters are available for this function.
#' For more information on how to implement the filtering parameters, refer to the parameter descriptions as well as examples in the vignettes.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_coding_ssm_status],
#' [GAMBLR.results::get_ssm_by_patients], [GAMBLR.results::get_ssm_by_sample], [GAMBLR.results::get_ssm_by_samples], [GAMBLR.results::get_ssm_by_region], [GAMBLR.results::get_ssm_by_regions]
#'
#' @param limit_cohort Supply this to restrict mutations to one or more cohorts in a vector.
#' @param exclude_cohort  Supply this to exclude mutations from one or more cohorts in a vector.
#' @param limit_pathology Supply this to restrict mutations to one pathology.
#' @param limit_samples This parameter is deprecated, use `these_sample_ids` instead.
#' @param these_samples_metadata Supply a metadata table to auto-subset the data to samples in that table before returning.
#' @param these_sample_id Optional, restrict the returned maf to a set of sample IDs (this parameter is replacing the `limit_samples` parameter).
#' @param force_unmatched_samples Optional argument for forcing unmatched samples, using [GAMBLR.results::get_ssm_by_samples].
#' @param projection Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#' @param this_seq_type The seq_type you want back, default is genome.
#' @param basic_columns Set to FALSE to override the default behavior of returning only the first 45 columns of MAF data.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters (matching columns in MAF).
#' @param from_flatfile Set to TRUE to obtain mutations from a local flatfile instead of the database. This can be more efficient and is currently the only option for users who do not have ICGC data access.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs)
#' @param groups Unix groups for the samples to be included. Default is both gambl and icgc_dart samples.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is TRUE.
#' @param engine Currently, only fread_maf (default) is supported.
#' @param verbose Controls the "verboseness" of this function (and internally called helpers).
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#'
#' @import dplyr tidyr RMariaDB DBI glue GAMBLR.helpers
#' @export
#'
#' @examples
#' #basic usage
#' maf_data = get_coding_ssm(limit_cohort = c("BL_ICGC"))
#'
#' maf_data = get_coding_ssm(limit_samples = "HTMCP-01-06-00485-01A-01D")
#'
get_coding_ssm = function(limit_cohort,
                          exclude_cohort,
                          limit_pathology,
                          limit_samples,
                          these_sample_ids = NULL,
                          these_samples_metadata = NULL,
                          force_unmatched_samples,
                          projection = "grch37",
                          this_seq_type = "genome",
                          basic_columns = TRUE,
                          maf_cols = NULL,
                          from_flatfile = TRUE,
                          augmented = TRUE,
                          min_read_support = 3,
                          groups = c("gambl", "icgc_dart"),
                          include_silent = TRUE,
                          engine = "fread_maf",
                          verbose = TRUE){
  
  if(!missing(limit_samples)){
    stop("limit_samples is deprecated, use `these_sample_ids` instead, or use a metadata table already subset to the samples of interest with `these_samples_metadata`")
  }

  remote_session = check_remote_configuration()
  
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  
  #get samples with the dedicated helper function
  metadata = id_ease(these_samples_metadata = these_samples_metadata,
                     these_sample_ids = these_sample_ids,
                     verbose = verbose,
                     this_seq_type = this_seq_type)
  
  #do all remaining filtering on the metadata then add the remaining sample_id to the query
  #unix groups
  metadata = metadata %>%
    dplyr::filter(unix_group %in% groups)

  #limit cohort
  if(!missing(limit_cohort)){
    metadata = metadata %>%
      dplyr::filter(cohort %in% limit_cohort)
  }
  #exclude cohort
  if(!missing(exclude_cohort)){
    metadata = metadata %>%
      dplyr::filter(!cohort %in% exclude_cohort)
  }
  #limit pathology
  if(!missing(limit_pathology)){
    metadata = metadata %>%
      dplyr::filter(pathology %in% limit_pathology)
  }
  
  #pull info for loading .CDS.maf
  sample_ids = pull(metadata, sample_id)

  #get file path for non-augmented maf
  if(from_flatfile && !augmented){
    seq_type = this_seq_type #glue needs to have this variable
    maf_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$ssm$template$cds$deblacklisted)
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), maf_path)
  }

  #get file path for augmented maf
  if(from_flatfile && augmented){
    seq_type = this_seq_type #glue needs to have this variable
    maf_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$ssm$template$cds$augmented)
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), maf_path)
  }

  #read file
  if(from_flatfile){
    if(verbose){
      message(paste("reading from:", full_maf_path)) 
    }

  #check for missingness
    if(!file.exists(full_maf_path)){
      print(paste("missing: ", full_maf_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    if(engine=='fread_maf'){
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
    }else{
      stop("Currently, only fread_maf is the only supported read engine...")
    }
    
    if(verbose){
      mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
      message(paste("mutations from", mutated_samples, "samples")) 
    }
    
  }else{
    #use db if not using flat file (mostly deprecated)...

    table_name = GAMBLR.helpers::check_config_value(config::get("results_tables")$ssm)
    db = GAMBLR.helpers::check_config_value(config::get("database_name"))
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

    muts = tbl(con, table_name) %>%
      dplyr::filter(Variant_Classification %in% coding_class) %>%
      as.data.frame()

    DBI::dbDisconnect(con)
  }

  #if augmented maf selected, drop variants with low read support (default is 3)
  if(augmented){
    muts = dplyr::filter(muts, t_alt_count >= min_read_support)
  }

  #filter maf on selected sample ids
  muts = muts %>%
    dplyr::filter(Tumor_Sample_Barcode %in% sample_ids)

  if(verbose){
    mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
    message(paste("after linking with metadata, we have mutations from", mutated_samples, "samples")) 
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
  return(muts)
}
