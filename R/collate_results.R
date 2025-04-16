#' @title Collate Results
#'
#' @description Bring together all derived sample-level results from many GAMBL pipelines.
#'
#' @details This function takes a data frame with sample IDs and seq types (in sample_id and seq_type columns) with the `these_samples_metadata` parameter and adds sample-level results from many of the available GAMBL pipelines.
#' Optional parameters are `join_with_full_metadata`. If `join_with_full_metadata` is set to TRUE, the function will return all metadata with `get_gambl_metadata`, allowing the user to extend the available information in a metadata table.
#' This function has also been designed so that it can get cached results, meaning that not all individual collate helper functions would have to be run to get results back.
#' To do so, run this function with `from_cache = TRUE` (default). In addition, it's also possible to regenerate the cached results, this is done by setting `write_to_file = TRUE`,
#' This operation auto defaults `from_cache = FALSE`. `case_set` is an optional parameter available for subsetting the return to an already defined set of cases.
#' If a dataframe is not provided, the function will default to all genome metadata returned with `get_gambl_metadata`. For more information on how to get the most out of this function,
#' refer to function examples, vignettes and parameter descriptions.
#'
#' @param sample_table A data frame with sample_id as the first column (deprecated, use these_samples_metadata instead).
#' @param write_to_file Boolean statement that outputs tsv file (/projects/nhl_meta_analysis_scratch/gambl/results_local/shared/gambl_{seq_type_filter}_results.tsv) if TRUE, default is FALSE.
#' @param join_with_full_metadata Join with all columns of metadata, default is FALSE.
#' @param these_samples_metadata Optional argument to use a user specified metadata df, must include seq_type column to specify desired seq types. If not provided, defaults to all genome samples with `get_gambl_metadata`. 
#' @param case_set Optional short name for a pre-defined set of cases.
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#' @param seq_type_filter Filtering criteria, default is genomes (deprecated. Include a `seq_type` column in these_samples_metadata to specify seq type).
#' @param from_cache Boolean variable for using cached results (/projects/nhl_meta_analysis_scratch/gambl/results_local/shared/gambl_{seq_type}_results.tsv), default is TRUE. If write_to_file is TRUE, this parameter auto-defaults to FALSE.
#'
#' @return A table keyed on biopsy_id that contains a bunch of per-sample results from GAMBL
#'
#' @import dplyr readr config glue GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' #get collated results for all capture samples, using cached results
#' capture_metadata = get_gambl_metadata(dna_seq_type_priority = "capture") %>% dplyr::filter(seq_type == "capture")
#' capture_collated_everything = collate_results(these_samples_metadata = capture_metadata,
#'                                               from_cache = TRUE,
#'                                               write_to_file = FALSE)
#'
#' #use an already subset metadata table for getting collated results (cached)
#' my_metadata = get_gambl_metadata()
#' fl_metadata = dplyr::filter(my_metadata, pathology == "FL")
#'
#' fl_collated = collate_results(these_samples_metadata = fl_metadata,
#'                               write_to_file = FALSE,
#'                               from_cache = TRUE)
#' }
#' #use an already subset metadata table for getting collated results (without using cached results)
#' my_metadata = get_gambl_metadata()
#' fl_metadata = dplyr::filter(my_metadata, pathology == "FL")
#'
#' fl_collated = collate_results(these_samples_metadata = fl_metadata,
#'                               write_to_file = FALSE,
#'                               from_cache = FALSE)
#' dplyr::select(fl_collated, 1:14) %>% head()
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
  
  if (!missing(sample_table)){
    stop("Argument sample_table is deprecated. Use these_samples_metadata instead")
  }
  if(!missing(seq_type_filter)){
    stop("Argument seq_type_filter is deprecated. ")
  }

  if(!missing(these_samples_metadata)){
    sample_table = these_samples_metadata
  } else {
    print("No sample table or metadata dataframe provided, defaulting to genome samples. Provide a dataframe that includes seq_type column to either the sample_table or these_samples_metadata argument to specify desired samples and seq type(s).")
    sample_table = get_gambl_metadata() %>%
      dplyr::filter(seq_type %in% c("genome")) %>% 
      dplyr::select(sample_id, patient_id, biopsy_id, seq_type)
  }

  if(!any(grepl("seq_type", names(sample_table)))){
    stop(paste0("Please specify desired seq type by including a seq_type column in these_samples_metadata dataframe"))
  }
  if(write_to_file){
    from_cache = FALSE #override default automatically for nonsense combination of options
  }

  seq_types = unique(sample_table$seq_type)
  seq_types_df = data.frame(seq_type = seq_types, cols = seq_types)

  #get paths to cached results, for from_cache = TRUE and for writing new cached results.
  output_file = GAMBLR.helpers::check_config_value(config::get("results_merged")$collated)
  output_base = GAMBLR.helpers::check_config_value(config::get("project_base"))
  output_file = paste0(output_base, output_file)
  output_file = lapply(seq_types, function(x) glue::glue(output_file, seq_type_filter = x))
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
    sample_table = do.call(bind_rows, lapply(
      output_file, 
      function(x) mutate(read_tsv(x), seq_type = case_when(
        grepl("genome", x) ~ "genome",
        grepl("capture", x) ~ "capture",
        grepl("mrna", x) ~ "mrna")))) %>% 
      right_join(., dplyr::select(sample_table, sample_id, patient_id, biopsy_id, seq_type))

  }else{
    message("Slow option: not using cached result. I suggest from_cache = TRUE whenever possible")
    #edit this function and add a new function to load any additional results into the main summary table
    seen_cols = c()
    original_cols = colnames(sample_table)
    sample_table_final = slice(sample_table, 0)
    for (seq in seq_types){
      sample_table_temp = collate_ssm_results(sample_table = dplyr::filter(sample_table, seq_type == seq), seq_type_filter = seq)
      sample_table_temp = collate_sv_results(sample_table = sample_table_temp, seq_type_filter = seq)
      sample_table_temp = collate_curated_sv_results(sample_table = sample_table_temp, seq_type_filter = seq)
      if(seq != "mrna") {
        sample_table_temp = collate_ashm_results(sample_table = sample_table_temp, seq_type_filter = seq)
        sample_table_temp = collate_nfkbiz_results(sample_table = sample_table_temp, seq_type_filter = seq)
      }
      #sample_table_temp = collate_csr_results(sample_table = sample_table_temp, seq_type_filter = seq)
      sample_table_temp = collate_ancestry(sample_table = sample_table_temp, seq_type_filter = seq)
      sample_table_temp = collate_sbs_results(sample_table = sample_table_temp, sbs_manipulation = sbs_manipulation, seq_type_filter = seq)
      if(seq != "mrna") {
        sample_table_temp = collate_qc_results(sample_table = sample_table_temp, seq_type_filter = seq)
      }
      #sample_table_temp = collate_pga(sample_table = sample_table_temp, this_seq_type = seq)
      sample_table_temp = collate_dlbclass(sample_table = (sample_table_temp))

      seq_types_df[seq_types_df$seq_type==seq,]$cols = list(colnames(sample_table_temp))

      # Handle incompatible column types
      if(length(seen_cols) != 0){
        common_cols = intersect(colnames(sample_table_temp), seen_cols)

        diff_class_cols = data.frame(
          Column = common_cols, 
          prev_class = sapply(sample_table_final[,common_cols], class),
          current_class = sapply(sample_table_temp[,common_cols], class)) %>%
          dplyr::filter(prev_class != current_class)

        for (col in pull(diff_class_cols, Column)) {
          row = diff_class_cols[diff_class_cols$Column == col,]
          
          if(row$prev_class %in% c("numeric","double","integer") & row$current_class == "character") {
            sample_table_temp[[col]] = as.numeric(sample_table_temp[[col]])
          }
          if(row$current_class %in% c("numeric","double","integer") & row$prev_class == "character") {
            sample_table_final[[col]] = as.numeric(sample_table_final[[col]])
          }
          if(row$prev_class == "logical" & row$current_class == "character") {
            sample_table_final[[col]] = as.character(sample_table_final[[col]])
          }
          if (row$current_class == "logical" & row$prev_class == "character") {
            sample_table_temp[[col]] = as.character(sample_table_temp[[col]])
          }
        }
      }

      sample_table_final = bind_rows(sample_table_final, sample_table_temp)
      seen_cols = setdiff(colnames(sample_table_final), original_cols)
    }
    sample_table = sample_table_final
  }
  if(write_to_file){
    #write results from "slow option" to new cached results file
    for (f in output_file) {
      if (grepl("genome", f)) {
        write_tsv(
          dplyr::select(dplyr::filter(sample_table, seq_type == "genome"), 
          all_of(seq_types_df[seq_types_df$seq_type=="genome",]$cols[[1]])), 
          file = f)
      }
      if (grepl("capture", f)) {
        write_tsv(
          dplyr::select(dplyr::filter(sample_table, seq_type == "capture"),
          all_of(seq_types_df[seq_types_df$seq_type=="capture",]$cols[[1]])), 
          file = f)
      }
      if (grepl("mrna", f)) {
        write_tsv(
          dplyr::select(dplyr::filter(sample_table, seq_type == "mrna"),
          all_of(seq_types_df[seq_types_df$seq_type=="mrna",]$cols[[1]])),
          file = f)
      }
    }
  }
  #convenience columns bringing together related information
  if(join_with_full_metadata){
    if (all(c("capture","genome") %in% seq_types)) {
      meta_data_genome = get_gambl_metadata(dna_seq_type_priority="genome")
      meta_data_capture = get_gambl_metadata(dna_seq_type_priority="capture")
      meta_data_capture = right_join(
        meta_data_capture, 
        anti_join(
          dplyr::select(meta_data_capture, sample_id, patient_id, biopsy_id, seq_type),
          dplyr::select(meta_data_genome, sample_id, patient_id, biopsy_id, seq_type)
          )
        )
      meta_data = bind_rows(meta_data_genome, meta_data_capture)
    } else {
      meta_data <- if(length(seq_types[!seq_types=="mrna"])==0){get_gambl_metadata()}else{get_gambl_metadata(dna_seq_type_priority=seq_types[!seq_types=="mrna"])}
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
