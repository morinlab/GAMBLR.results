#' @title Get SSM By Samples.
#'
#' @description Get MAF-format data frame for more than one sample and combine them together.
#'
#' @details This function internally runs [GAMBLR.results::get_ssm_by_sample].
#' The user can either give the function a vector of sample IDs of interest with `these_sample_ids`,
#' or use a metadata table (`these_samples_metadata`), already subset to the sample IDs of interest.
#' In most situations, this should never need to be run with subset_from_merge = TRUE.
#' Instead use one of [GAMBLR.results::get_coding_ssm] or [GAMBLR.results::get_ssm_by_region].
#' See [GAMBLR.results::get_ssm_by_sample] for more information.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_coding_ssm],
#' [GAMBLR.utils::get_coding_ssm_status], [GAMBLR.results::get_ssm_by_patients], [GAMBLR.results::get_ssm_by_sample], [GAMBLR.results::get_ssm_by_region], [GAMBLR.results::get_ssm_by_regions]
#'
#' @param these_sample_ids A vector of sample_id that you want results for.
#' @param these_samples_metadata Optional metadata table. If provided, the function will return SSM calls for the sample IDs in the provided metadata table.
#' @param tool_name Only supports slms-3 currently.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param seq_type  The seq type you want results for. Default is "genome".
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatibility with additional variant calling parameters/versions.
#' @param these_genes A vector of genes to subset ssm to.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 45 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters.
#' @param subset_from_merge Instead of merging individual MAFs, the data will be subset from a pre-merged MAF of samples with the specified seq_type.
#' @param engine Specify one of readr or fread_maf (default) to change how the large files are loaded prior to subsetting. You may have better performance with one or the other but for me fread_maf is faster and uses a lot less RAM.
#'
#' @return A data frame in MAF format.
#'
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' #examples using the these_sample_ids parameter.
#' sample_ssms = get_ssm_by_samples(these_sample_ids = c("HTMCP-01-06-00485-01A-01D",
#'                                                       "14-35472_tumorA",
#'                                                       "14-35472_tumorB"))
#'
#' hg38_ssms = get_ssm_by_samples(projection="hg38",
#'                                these_sample_ids = c("HTMCP-01-06-00485-01A-01D",
#'                                                     "14-35472_tumorA",
#'                                                     "14-35472_tumorB"))
#'
#' readr_sample_ssms = get_ssm_by_samples(subset_from_merge = TRUE,
#'                                        engine = "readr",
#'                                        these_sample_ids = c("HTMCP-01-06-00485-01A-01D",
#'                                                             "14-35472_tumorA",
#'                                                             "14-35472_tumorB"))
#'
#' slow_sample_ssms = get_ssm_by_samples(subset_from_merge = TRUE,
#'                                       these_sample_ids = c("HTMCP-01-06-00485-01A-01D",
#'                                                            "14-35472_tumorA",
#'                                                            "14-35472_tumorB"))
#'
#' #example using a metadata table subset to sample IDs of interest.
#' my_metadata = get_gambl_metadata(seq_type_filter = "genome")
#' my_metadata = dplyr::filter(my_metadata, pathology == "FL")
#'
#' sample_ssms = get_ssm_by_samples(these_samples_metadata = my_metadata)
#'
get_ssm_by_samples = function(these_sample_ids,
                              these_samples_metadata,
                              tool_name = "slms-3",
                              projection = "grch37",
                              seq_type = "genome",
                              flavour = "clustered",
                              these_genes,
                              min_read_support = 3,
                              basic_columns = TRUE,
                              maf_cols = NULL,
                              subset_from_merge = FALSE,
                              augmented = TRUE,
                              engine = 'fread_maf'){

  remote_session = check_remote_configuration(auto_connect = TRUE)
  if(!subset_from_merge){
    message("WARNING: on-the-fly merges can be extremely slow and consume a lot of memory if many samples are involved. Use at your own risk. ")
  }
  to_exclude = get_excluded_samples(tool_name)

  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = seq_type) %>%
      dplyr::filter(sample_id %in% these_sample_ids) %>%
      dplyr::filter(!sample_id %in% to_exclude)
  }else{
    if(missing(these_sample_ids)){
      #assume the user just wants the data for all the sample ids in this data frame
      these_sample_ids = pull(these_samples_metadata,sample_id)
    }
    else{
      these_samples_metadata = these_samples_metadata %>%
        dplyr::filter(sample_id %in% these_sample_ids) %>%
        dplyr::filter(!sample_id %in% to_exclude)
    }
  }
  #ensure we only have sample_id that are in the remaining metadata (no excluded/unavailable samples)
  these_sample_ids = these_sample_ids[which(these_sample_ids %in% these_samples_metadata$sample_id)]
  maf_column_types = "ccccciiccccccccccccccccccccccnccccccccciiiiii" #for the first 45 standard columns
  if(flavour=="legacy"){
    warning("I lied. Access to the old variant calls is not currently supported in this function")
    # TODO: implement loading of the old merged MAF under icgc_dart... vcf2maf-1.2 ..level_3 as per the other from_flatfile functions
    return()

  }else if(flavour=="clustered"){
    if(subset_from_merge && !augmented){
      maf_template = check_config_value(config::get("results_flatfiles")$ssm$template$merged$deblacklisted)
      maf_path = glue::glue(maf_template)
      full_maf_path =  paste0(check_config_value(config::get("project_base")), maf_path)
      message(paste("using existing merge:", full_maf_path))

      #check for missingness
      if(!file.exists(full_maf_path)){
        print(paste("missing: ", full_maf_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
        }

      if(engine=="fread_maf"){
        if(basic_columns){
          maf_df_merge = fread_maf(full_maf_path,select_cols = c(1:45)) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }else{
          maf_df_merge = fread_maf(full_maf_path) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }
      }else if(engine=="readr"){
        if(basic_columns){
          maf_df_merge = suppressMessages(read_tsv(full_maf_path,col_select = c(1:45),num_threads=12,col_types = maf_column_types,lazy = TRUE)) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }else{
          maf_df_merge = fread_maf(full_maf_path) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }
      }else{
        stop("specify one of readr or fread_maf as the file-reading engine")
      }

      if(!is.null(maf_cols) && !basic_columns){
        maf_df_merge = dplyr::select(maf_df_merge, all_of(maf_cols))
      }
    }

    if(subset_from_merge && augmented){
      maf_template = check_config_value(config::get("results_flatfiles")$ssm$template$merged$augmented)
      maf_path = glue::glue(maf_template)
      full_maf_path =  paste0(check_config_value(config::get("project_base")), maf_path)
      message(paste("using existing merge:", full_maf_path))

      #check for missingness
      if(!file.exists(full_maf_path)){
        print(paste("missing: ", full_maf_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
      }

      #maf_df_merge = read_tsv(full_maf_path) %>%
      #  dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
      #  dplyr::filter(t_alt_count >= min_read_support)
      if(basic_columns){
        maf_df_merge = fread_maf(full_maf_path,select_cols = c(1:45)) %>%
          dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
          dplyr::filter(t_alt_count >= min_read_support)
      }else{
        maf_df_merge = fread_maf(full_maf_path) %>%
          dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
          dplyr::filter(t_alt_count >= min_read_support)
      }
      #subset maf to only include first 43 columns (default)
      if(basic_columns){maf_df_merge = dplyr::select(maf_df_merge, c(1:45))}
      #subset maf to a specific set of columns (defined in maf_cols)
      if(!is.null(maf_cols) && !basic_columns){maf_df_merge = dplyr::select(maf_df_merge, all_of(maf_cols))}
    }

    if(!subset_from_merge){
      if(remote_session){
        maf_df_list = list()
        for(this_sample in these_sample_ids){
          maf_df = get_ssm_by_sample(
            this_sample_id = this_sample,
            these_samples_metadata = these_samples_metadata,
            tool_name = tool_name,
            projection = projection,
            augmented = augmented,
            flavour = flavour,
            min_read_support = min_read_support,
            basic_columns = basic_columns,
            maf_cols = maf_cols,
            verbose = FALSE)
          maf_df_list[[this_sample]]=maf_df
        }
      }else{
        maf_df_list = parallel::mclapply(these_sample_ids,function(x){get_ssm_by_sample(
        this_sample_id=x,
        these_samples_metadata = these_samples_metadata,
        tool_name = tool_name,
        projection = projection,
        augmented = augmented,
        flavour = flavour,
        min_read_support = min_read_support,
        basic_columns = basic_columns,
        maf_cols = maf_cols,
        verbose = FALSE
        )},mc.cores = 12)
      }
      maf_df_merge = bind_rows(maf_df_list)
    }
  }

  if(!missing(these_genes)){
    maf_df_merge = maf_df_merge %>%
      dplyr::filter(Hugo_Symbol %in% these_genes)
  }
  return(maf_df_merge)
}
