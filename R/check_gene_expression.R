
#' @title Check Gene Expression.
#'
#' @description This function determines which samples have expression data available in the merge and drop redundant data while consistently prioritizing by protocol and nucleic acid source.
#' 
#' 
#' @param these_samples_metadata The data frame with sample metadata. Usually output of the get_gambl_metadata().
#' @param verbose Set to TRUE mainly for debugging
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
check_gene_expression = function(verbose=F){
  # We start with the minimal metadata for all available samples in the merge.
  # This file allows you to load the required information without loading the entire tidy expression file
  # TODO: put this in the config and automate its generation when the main files are created so it remains up-to-date
  onegene = paste0("/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/DESeq2-0.0_salmon-1.0/mrna--gambl-icgc-all/","vst-just_onegene.tsv")
  
  expression_data_rows = suppressMessages(read_tsv(onegene,col_select = c(mrna_sample_id,patient_id,biopsy_id,protocol,ffpe_or_frozen),progress = F)) %>% 
    mutate(Available=TRUE) %>%
    mutate(seq_type="mrna")
  these_samples_metadata = get_gambl_metadata(seq_type_filter = c("genome","capture"))  %>%
    dplyr::select(sample_id, 
                  patient_id, biopsy_id, seq_type)
  
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
  
  
  # collapse duplicates if any
  if(any(expression_data_rows$multi_exp == 1) ){
    if(verbose){
      rows_total = filter(expression_data_rows,multi_exp==1) %>% nrow()
      unique_redundant = filter(expression_data_rows,multi_exp==1) %>% pull(biopsy_id) %>% unique()
      nid = length(unique_redundant)
      print(paste(rows_total, "rows and ", nid, "unique biopsies are redundant, dropping redundancy"))
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
  expression_data_rows = filter(expression_data_rows,Available==TRUE) %>% 
    select(-Available) %>%
    mutate(sample_id = mrna_sample_id) %>% 
    unique()

  #join to other seq_types where possible
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
