#' @title Get SSM By Samples.
#'
#' @description Get the genome-wide set of mutations for one or more sample including coding and non-coding mutations.
#'
#' @details 
#' The user can specify a metadata table (`these_samples_metadata`), subset to the sample IDs of interest.
#' In most situations, this should never need to be run with subset_from_merge = TRUE, which is very inefficient.
#' This function does not scale well to many samples. In most cases, users will actually need either [GAMBLR.results::get_coding_ssm] or [GAMBLR.results::get_ssm_by_region].
#' See [GAMBLR.results::get_ssm_by_sample] for more information.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_coding_ssm],
#' [GAMBLR.results::get_ssm_by_patients], [GAMBLR.results::get_ssm_by_regions]
#'
#' @param these_samples_metadata Optional metadata table.
#' If provided, it will return SSM calls for the samples in the metadata table.
#' @param tool_name Only supports slms-3 currently.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 45 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters.
#' @param subset_from_merge Instead of merging individual MAFs, the data will be subset from a pre-merged MAF of samples with the specified this_seq_type.
#' @param engine Specify one of readr or fread_maf (default) to change how the large files are loaded prior to subsetting. You may have better performance with one or the other but for me fread_maf is faster and uses a lot less RAM.
#' @param this_seq_type Deprecated. Inferred from these_samples_metadata
#' @param these_sample_ids Deprecated. Inferred from these_samples_metadata
#'
#' @return A data frame in MAF format.
#'
#' @import dplyr readr tidyr glue GAMBLR.helpers
#' @export
#'
#' @examples
#' 
#' my_meta = get_gambl_metadata() %>% 
#'                        dplyr::filter(sample_id %in% c("HTMCP-01-06-00485-01A-01D",
#'                                                "14-35472_tumorA",
#'                                                "14-35472_tumorB"))
#' sample_ssms = get_ssm_by_samples(these_samples_metadata = my_meta)
#'
#' hg38_ssms = get_ssm_by_samples(projection="hg38",
#'                                these_samples_metadata = my_meta)
#' 
#' dplyr::group_by(hg38_ssms,Tumor_Sample_Barcode) %>% count()
#' hg38_ssms_no_aug = get_ssm_by_samples(projection="hg38",
#'                                these_samples_metadata = my_meta,augmented= FALSE)
#' 
#' dplyr::group_by(hg38_ssms_no_aug,Tumor_Sample_Barcode) %>% count()
#' \dontrun{
#' my_metadata = dplyr::filter(my_metadata, pathology == "FL")
#'
#' sample_ssms = get_ssm_by_samples(these_samples_metadata = my_metadata)
#' }
get_ssm_by_samples = function(these_samples_metadata,
                              tool_name = "slms-3",
                              projection = "grch37",
                              flavour = "clustered",
                              these_genes,
                              min_read_support = 3,
                              basic_columns = TRUE,
                              maf_cols = NULL,
                              subset_from_merge = FALSE,
                              augmented = TRUE,
                              engine = 'fread_maf',
                              these_sample_ids,
                              this_seq_type){
  
  if(!missing(this_seq_type) | !missing(these_sample_ids)){
    stop("this_seq_type and these_sample_ids are deprecated. Use these_samples_metadata instead")
  }
  remote_session = check_remote_configuration(auto_connect = TRUE)

  to_exclude = get_excluded_samples(tool_name)

  if(missing(these_samples_metadata)){
    warning(
        "CAUTION! The option `these_sample_ids` may behave unexpectedly if a mix of seq_types is expected."
    )
    warning(
        "Please use the `these_samples_metadata` option instead of `these_sample_ids`."
    )
    these_samples_metadata = get_gambl_metadata() %>%
      dplyr::filter(sample_id %in% these_sample_ids,seq_type!="mrna") %>%
      dplyr::filter(!sample_id %in% to_exclude)
      n_samp = nrow(these_samples_metadata)
      n_orig = length(these_sample_ids)
      if(n_orig > n_samp){
        lost = n_orig - n_samp
        print(paste(lost,"of",n_orig,"samples have been dropped!"))
      }
  }else{
   #drop unsupported seq_type
    these_samples_metadata = dplyr::filter(these_samples_metadata,seq_type!="mrna")
    if(missing(these_sample_ids)){
      #assume the user just wants the data for all the sample ids in this data frame
      seq_type_sample_ids = list()
      for(a_seq_type in unique(these_samples_metadata$seq_type)){
        seq_type_sample_ids[[a_seq_type]]=dplyr::filter(these_samples_metadata,seq_type==a_seq_type) %>% pull(sample_id)
      } 
    }
    else{
      these_samples_metadata = these_samples_metadata %>%
        dplyr::filter(sample_id %in% these_sample_ids) %>%
        dplyr::filter(!sample_id %in% to_exclude)
    }
    if(length(unique(these_samples_metadata$seq_type))>1){
      message("metadata provided contained more than one seq_type.")
      print("The value of this_seq_type was ignored and all samples in these_samples_metadata and these_sample_ids (if provided) were included")
    }
  }
  #ensure we only have sample_id that are in the remaining metadata (no excluded/unavailable samples)
  
  maf_column_types = "ccccciiccccccccccccccccccccccnccccccccciiiiii" #for the first 45 standard columns
  if(flavour=="legacy"){
    warning("I lied. Access to the old variant calls is not currently supported in this function")
    # TODO: implement loading of the old merged MAF under icgc_dart... vcf2maf-1.2 ..level_3 as per the other from_flatfile functions
    return()

  }else if(flavour=="clustered"){
    if(subset_from_merge && !augmented){
      if(length(unique(these_samples_metadata$seq_type))>1){
        print("more than one seq_type provided")
        print("This function needs to be updated to handle >1 seq_type")
        print("For a workaround, you can run individually for each desired seq_type")
        stop()
      }
      seq_type = these_samples_metadata$seq_type[1] #needed for glue
      maf_template = check_config_and_value("results_flatfiles$ssm$template$merged$deblacklisted")
      maf_path = glue::glue(maf_template)
      full_maf_path =  paste0(check_config_and_value("project_base"), maf_path)
      message(paste("using existing merge:", full_maf_path))
      #if(!file.exists(full_maf_path)){
      #  full_maf_path = paste0(full_maf_path,".bgz")  
      #}
      #check for missingness
      if(!file.exists(full_maf_path)){
        print(paste("missing: ", full_maf_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
      }

      
      if(engine=="fread_maf"){
        if(basic_columns){
          maf_df_merge = suppressMessages(fread_maf(full_maf_path,select_cols = c(1:45))) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }else{
          maf_df_merge = suppressMessages(fread_maf(full_maf_path)) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }
      }else if(engine=="readr"){
        if(basic_columns){
          maf_df_merge = suppressMessages(
            read_tsv(full_maf_path,
            col_select = c(1:45),
            num_threads=12,col_types = maf_column_types,lazy = TRUE)) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }else{
          maf_df_merge = suppressMessages(fread_maf(full_maf_path)) %>%
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
      if(length(unique(these_samples_metadata$seq_type))>1){
        print("more than one seq_type provided")
        print("This function needs to be updated to handle >1 seq_type")
        print("For a workaround, you can run individually for each desired seq_type")
        stop()
      }
      
      seq_type = these_samples_metadata$seq_type[1] #needed for glue
      maf_template = check_config_and_value("results_flatfiles$ssm$template$merged$augmented")
      maf_path = glue::glue(maf_template)
      full_maf_path =  paste0(check_config_and_value("project_base"), maf_path)
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
        maf_df_merge = suppressMessages(fread_maf(full_maf_path,select_cols = c(1:45))) %>%
          dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
          dplyr::filter(t_alt_count >= min_read_support)
      }else{
        maf_df_merge = suppressMessages(fread_maf(full_maf_path)) %>%
          dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
          dplyr::filter(t_alt_count >= min_read_support)
      }
      #subset maf to only include first 43 columns (default)
      if(basic_columns){maf_df_merge = dplyr::select(maf_df_merge, c(1:45))}
      #subset maf to a specific set of columns (defined in maf_cols)
      if(!is.null(maf_cols) && !basic_columns){maf_df_merge = dplyr::select(maf_df_merge, all_of(maf_cols))}
    }

    if(!subset_from_merge){
      maf_df_list = list()
      for(a_seq_type in names(seq_type_sample_ids)){
          for(this_sample in seq_type_sample_ids[[a_seq_type]]){
            maf_df = get_ssm_by_sample(
              these_samples_metadata = dplyr::filter(these_samples_metadata,
                                                     sample_id==this_sample,
                                                     seq_type==a_seq_type),
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
        }
      
      maf_df_merge <- do.call(bind_genomic_data, maf_df_list)
    }
  }

  return(maf_df_merge)
}
