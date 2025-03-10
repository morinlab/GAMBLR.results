#' @title Collate Results
#'
#' @description Bring together all derived sample-level results from many GAMBL pipelines.
#'
#' @details This function takes a data frame with sample IDs (in the first column) with the `sample_table` parameter and adds sample-level results from many of the available GAMBL pipelines.
#' Optional parameters are `these_samples_metadata` and `join_with_full_metadata`. If `join_with_full_metadata` is set to TRUE, the function can either work with an already subset metadata
#' table (`these_sampels_metadata`), or, if not provided, the function will default to all metadata returned with `get_gambl_metadata`, allowing the user to extend the available information in a metadata table.
#' This function has also been designed so that it can get cached results, meaning that not all individual collate helper functions would have to be run to get results back.
#' To do so, run this function with `from_cache = TRUE` (default). In addition, it's also possible to regenerate the cached results, this is done by setting `write_to_file = TRUE`,
#' This operation auto defaults `from_cache = FALSE`. `case_set` is an optional parameter available for subsetting the return to an already defined set of cases.
#' Lastly, `seq_type_filter` lets the user control what seq type results will be returned for. Default is "genome". For more information on how to get the most out of this function,
#' refer to function examples, vignettes and parameter descriptions.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param write_to_file Boolean statement that outputs tsv file (/projects/nhl_meta_analysis_scratch/gambl/results_local/shared/gambl_{seq_type_filter}_results.tsv) if TRUE, default is FALSE.
#' @param join_with_full_metadata Join with all columns of metadata, default is FALSE.
#' @param these_samples_metadata Optional argument to use a user specified metadata df, overwrites get_gambl_metadata in join_with_full_metadata.
#' @param case_set Optional short name for a pre-defined set of cases.
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param from_cache Boolean variable for using cached results (/projects/nhl_meta_analysis_scratch/gambl/results_local/shared/gambl_{seq_type_filter}_results.tsv), default is TRUE. If write_to_file is TRUE, this parameter auto-defaults to FALSE.
#'
#' @return A table keyed on biopsy_id that contains a bunch of per-sample results from GAMBL
#'
#' @import dplyr readr config glue GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' #get collated results for all capture samples, using cached results
#' capture_collated_everything = collate_results(seq_type_filter = "capture",
#'                                               from_cache = TRUE,
#'                                               write_to_file = FALSE)
#'
#' #use an already subset metadata table for getting collated results (cached)
#' my_metadata = get_gambl_metadata()
#' fl_metadata = dplyr::filter(my_metadata, pathology == "FL")
#'
#' fl_collated = collate_results(seq_type_filter = "genome",
#'                               join_with_full_metadata = TRUE,
#'                               these_samples_metadata = fl_metadata,
#'                               write_to_file = FALSE,
#'                               from_cache = TRUE)
#'
#' #get collated results for all genome samples and join with full metadata
#' everything_collated = collate_results(seq_type_filter = "genome",
#'                                       from_cache = TRUE,
#'                                       join_with_full_metadata = TRUE)
#'
#' #another example demonstrating correct usage of the sample_table parameter.
#' fl_samples = dplyr::select(fl_metadata, sample_id, patient_id, biopsy_id)
#'
#' fl_collated = collate_results(sample_table = fl_samples,
#'                               seq_type_filter = "genome",
#'                               from_cache = TRUE)
#' }
collate_results = function(sample_table,
                           write_to_file = FALSE,
                           join_with_full_metadata = FALSE,
                           these_samples_metadata,
                           case_set,
                           sbs_manipulation = "",
                           seq_type_filter = "genome",
                           from_cache = TRUE){

  # important: if you are collating results from anything but WGS (e.g RNA-seq libraries) be sure to use biopsy ID as the key in your join
  # the sample_id should probably not even be in this file if we want this to be biopsy-centric
  
  if(!missing(these_samples_metadata)){
    if(!missing(sample_table)){
      print("Prioritizing these_samples_metadata as sample input")
      sample_table = these_samples_metadata
    }
  } else if (missing(sample_table)){
    print("Defaulting to genome and capture metadata. Pass a dataframe that includes seq_type column as either the sample_table or these_samples_metadata argument to specify desired seq type")
    sample_table = get_gambl_metadata() %>%
      dplyr::filter(seq_type %in% c("genome","capture")) %>% 
      dplyr::select(sample_id, patient_id, biopsy_id, seq_type)
  }
  
  dep_msg=""
  if(!missing(seq_type_filter)){
    dep_msg = "Argument seq_type_filter is deprecated. "
  }
  if(!any(grepl("seq_type", names(sample_table)))){
    print(paste0(dep_msg, "Please provide a dataframe that includes seq_type column to either the sample_table or these_samples_metadata argument to specify desired seq type"))
  }
  if(write_to_file){
    from_cache = FALSE #override default automatically for nonsense combination of options
  }

  #get paths to cached results, for from_cache = TRUE and for writing new cached results.
  output_file = GAMBLR.helpers::check_config_value(config::get("results_merged")$collated)
  output_base = GAMBLR.helpers::check_config_value(config::get("project_base"))
  output_file = paste0(output_base, output_file)
  output_file = lapply(unique(sample_table$seq_type), function(x) glue::glue(output_file, seq_type_filter=x))
  print(output_file)
  if(from_cache){
    #check for missingness
    missing_cache = sapply(output_file, function(x) !file.exists(x))
    if(any(missing_cache)){
      print(paste("missing: ", output_file[missing_cache]))
      message("Cannot find file(s) locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    #read cached results
    sample_table = do.call(bind_rows, lapply(output_file, function(x) read_tsv(x))) %>% 
      left_join(., select(sample_table, sample_id, patient_id, biopsy_id, seq_type)) %>% # retain seq_type column
      dplyr::filter(sample_id %in% sample_table$sample_id)

  }else{
    message("Slow option: not using cached result. I suggest from_cache = TRUE whenever possible")
    #edit this function and add a new function to load any additional results into the main summary table
    sample_table = do.call(bind_rows, lapply(unique(sample_table$seq_type), function(x) collate_ssm_results(sample_table = filter(sample_table, seq_type==x), seq_type_filter = x)))
    sample_table = collate_sv_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_curated_sv_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_ashm_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_nfkbiz_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    #sample_table = collate_csr_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_ancestry(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_sbs_results(sample_table = sample_table, sbs_manipulation = sbs_manipulation, seq_type_filter = seq_type_filter)
    sample_table = collate_qc_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    #sample_table <- collate_pga(
    #    these_samples_metadata = sample_table,
    #    this_seq_type = seq_type_filter
    #)
    sample_table <- collate_dlbclass(
        sample_table
    )
  }
  if(write_to_file){
    #write results from "slow option" to new cached results file
    write_tsv(sample_table, file = output_file)
  }
  #convenience columns bringing together related information
  if(join_with_full_metadata){
    if(!missing(these_samples_metadata)){
      meta_data = these_samples_metadata
    }else{
      meta_data = get_gambl_metadata(seq_type_filter = seq_type_filter)
    }

    full_table = left_join(meta_data, sample_table)

    #check for missing columns, add any missing columns and fill with NA
    col_check = c("ashm_MYC", "manta_MYC_sv", "ICGC_MYC_sv", "myc_ba", "ashm_BCL2", "manta_BCL2_sv", "ICGC_BCL2_sv", "bcl2_ba") #create a vector of the columns to check for
    missing = setdiff(col_check, names(full_table)) #return the missing columns
    full_table[missing] = NA #add missing columns to full_table and fill such columns with NAs

    full_table = full_table %>%
      mutate("MYC_SV_any" = case_when(ashm_MYC > 3 ~ "POS", manta_MYC_sv == "POS" ~ "POS", ICGC_MYC_sv == "POS" ~ "POS", myc_ba == "POS" ~ "POS", TRUE ~ "NEG"))

    full_table = full_table %>%
      mutate("BCL2_SV_any" = case_when(ashm_BCL2 > 3 ~ "POS", manta_BCL2_sv == "POS" ~ "POS", ICGC_BCL2_sv == "POS" ~ "POS", bcl2_ba == "POS" ~ "POS", TRUE ~ "NEG"))

    full_table = full_table %>%
      mutate("DoubleHitBCL2" = ifelse(BCL2_SV_any == "POS" & MYC_SV_any == "POS", "Yes", "No"))
    return(full_table)
  }
  return(sample_table)
}
