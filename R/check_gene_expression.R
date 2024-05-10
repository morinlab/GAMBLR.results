
#' @title Check Gene Expression.
#'
#' @description This function determines which samples have expression data available in the merge and drop redundant data while consistently prioritizing by protocol and nucleic acid source.
#' 
#' @param show_linkages Set to TRUE to link every row to an available capture and genome sample 
#' using get_gambl_metadata to prioritize each to at most one per biopsy_id
#' @param verbose Set to TRUE mainly for debugging
#' @param ... Optional parameters to pass along to get_gambl_metadata (only used if show_linkages = TRUE)
#' 
#' @return A data frame with a row for each non-redundant RNA-seq result and the following columns:
#' 
#' \describe{
#'   \item{mrna_sample_id}{The unique sample_id value that will match a single row from the GAMBL metadata where seq_type is mrna. }
#'   \item{biopsy_id}{The unique identifier for the source of nucleic acids.}
#'   \item{sample_id}{Identical to mrna_sample_id}
#'   \item{capture_sample_id}{When this biopsy has capture/exome data in the GAMBL metadata, the value will be the sample_id for that data. NA otherwise.}
#'   \item{genome_sample_id}{When this biopsy has genome data in the GAMBL metadata, the value will be the sample_id for that data. NA otherwise.}
#'   \item{patient_id}{The anonymized unique identifier for this patient. For BC samples, this will be Res ID.}
#'   \item{seq_type}{The assay type used to produce this data (will always be "mrna" in this case)}
#'   \item{protocol}{Specifies the RNA-seq library construction protocol.}
#'   \item{ffpe_or_frozen}{Specifies the way the source of nucleic acids was preserved. Either FFPE or frozen.}}
#' 
check_gene_expression = function(verbose=F, show_linkages=F, ...){
  # We start with the minimal metadata that corresponds to the tidy expression file (sample_metadata.tsv)
  metadata_file = check_config_value(config::get("results_merged")$tidy_expression_metadata)
  metadata_file = paste0(check_config_value(config::get("project_base")),metadata_file)
  expression_data_rows = suppressMessages(read_tsv(metadata_file)) %>% 
                                            dplyr::select(sample_id,tissue_status,seq_type,patient_id,biopsy_id,protocol,ffpe_or_frozen,cohort) %>%
                                            dplyr::rename(c("mrna_sample_id"="sample_id"))
  

  
  collapse_duplicates = TRUE

  if(verbose){
    nrows = nrow(expression_data_rows)
    print(paste("starting with",nrows,"rows"))
  }

  
  mark_duplicates_biopsy = function(ew){
    mutate(expression_data_rows, sample_seqType = paste(biopsy_id, seq_type)) %>% 
      filter(!is.na(mrna_sample_id)) %>% 
      with( split(mrna_sample_id, sample_seqType) ) %>% 
      .[lengths(.) > 1] %>% 
      lapply(unique) %>% 
      .[lengths(.) > 1] %>% 
      names %>% 
      sub(" .+", "", .) %>% 
      { mutate(expression_data_rows, multi_exp = ifelse(biopsy_id %in% ., 1, 0)) }
  }
  
  expression_data_rows = mark_duplicates_biopsy(expression_data_rows)
  
  if(verbose){
    now_row = nrow(expression_data_rows)
    print(paste(now_row,"remain after collapsing on biopsy_id"))
  }
  
  # collapse duplicates if any
  if(any(expression_data_rows$multi_exp == 1) ){
    if(verbose){
      rows_total = filter(expression_data_rows,multi_exp==1) %>% nrow()
      unique_redundant = filter(expression_data_rows,multi_exp==1) %>% pull(biopsy_id) %>% unique()
      nid = length(unique_redundant)
      print(paste(rows_total, "LINE 74: rows corresponding to ", nid, "unique biopsies are redundant, dropping redundancy"))
      print(head(unique_redundant))
    }
    original = expression_data_rows %>% filter(!is.na(mrna_sample_id))
    expression_data_rows = group_by(expression_data_rows, patient_id, biopsy_id) %>% 
      slice_max(str_detect(protocol, "Ribo"), n=1, with_ties = TRUE) %>% 
      slice_max(ffpe_or_frozen == "frozen", n=1, with_ties = TRUE) %>% 
      ungroup()
    
    
    expression_data_rows = mark_duplicates_biopsy(expression_data_rows)
    
    if(verbose){
      now_row = nrow(expression_data_rows)
      print(paste(now_row,"remain after collapsing on biopsy_id"))
    }

  }else{
    lost = data.frame()
  }
  expression_data_rows = 
    mutate(expression_data_rows,sample_id = mrna_sample_id) %>% 
    unique()
  if(!show_linkages){
    return(select(expression_data_rows,-multi_exp))
  }
  #Optionally join to other seq_types where possible in case a user wants to see which samples have both RNA and a source of mutations
  these_samples_metadata = get_gambl_metadata(...)  %>%
    dplyr:: filter(seq_type %in% c("genome","capture")) %>%
    dplyr::select(sample_id, 
                  patient_id, biopsy_id, seq_type)
  capture_meta = dplyr::filter(these_samples_metadata,
                               seq_type=="capture",
                               patient_id %in% expression_data_rows$patient_id) %>%
    select(sample_id,patient_id,biopsy_id) %>%
    rename(c("capture_sample_id"="sample_id"))
  expression_data_rows = left_join(expression_data_rows,capture_meta,by=c("patient_id","biopsy_id"))
  
  genome_meta =  dplyr::filter(these_samples_metadata,
                               seq_type=="genome",
                               patient_id %in% expression_data_rows$patient_id) %>%
    select(sample_id,patient_id,biopsy_id) %>%
    rename(c("genome_sample_id"="sample_id"))
  expression_data_rows = left_join(expression_data_rows,genome_meta, by=c("patient_id","biopsy_id"))
  
  return(select(expression_data_rows,-multi_exp))
}
