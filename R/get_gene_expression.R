
#' @title Check Gene Expression.
#'
#' @description This function determines which samples have expression data available in the merge and drop redundant data while consistently prioritizing by protocol and nucleic acid source.
#' 
#' 
#' @param these_samples_metadata The data frame with sample metadata. Usually output of the get_gambl_metadata().
#' @param verbose Set to TRUE mainly for debugging
#' 
#' 
check_gene_expression = function(these_samples_metadata,
                                 verbose=F){
  onegene = paste0("/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/DESeq2-0.0_salmon-1.0/mrna--gambl-icgc-all/","vst-just_onegene.tsv")
  expression_data_rows = suppressMessages(read_tsv(onegene,col_select = c(mrna_sample_id,patient_id,biopsy_id,protocol,ffpe_or_frozen),progress = F)) %>% 
    mutate(Available=TRUE) %>%
    mutate(seq_type="mrna")
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = c("genome","capture")) 
  }
  
  collapse_duplicates = TRUE
  if(verbose){
    nrows = nrow(expression_data_rows)
    print(paste("starting with",nrows,"rows"))
  }
  these_samples_metadata = dplyr::select(these_samples_metadata, sample_id, 
                                         patient_id, biopsy_id, seq_type)
  expression_wider = expression_data_rows
  
  mark_duplicates_biopsy = function(ew){
    mutate(expression_wider, sample_seqType = paste(biopsy_id, seq_type)) %>% 
      filter(!is.na(mrna_sample_id)) %>% 
      with( split(mrna_sample_id, sample_seqType) ) %>% 
      .[lengths(.) > 1] %>% 
      lapply(unique) %>% 
      .[lengths(.) > 1] %>% 
      names %>% 
      sub(" .+", "", .) %>% 
      { mutate(expression_wider, multi_exp = ifelse(biopsy_id %in% ., 1, 0)) }
  }
  
  expression_wider = mark_duplicates_biopsy(expression_wider)
  
  
  # collapse duplicates if any
  if( collapse_duplicates & any(expression_wider$multi_exp == 1) ){
    if(verbose){
      rows_total = filter(expression_wider,multi_exp==1) %>% nrow()
      unique_redundant = filter(expression_wider,multi_exp==1) %>% pull(biopsy_id) %>% unique()
      nid = length(unique_redundant)
      print(paste(rows_total, "rows and ", nid, "unique biopsies are redundant, dropping redundancy"))
      print(head(unique_redundant))
    }
    original = expression_wider %>% filter(!is.na(mrna_sample_id))
    expression_wider = group_by(expression_wider, patient_id, biopsy_id) %>% 
      slice_max(str_detect(protocol, "Ribo"), n=1, with_ties = TRUE) %>% 
      slice_max(ffpe_or_frozen == "frozen", n=1, with_ties = TRUE) %>% 
      ungroup()
    
    
    expression_wider = mark_duplicates_biopsy(expression_wider)
    
    if(verbose){
      now_row = nrow(expression_wider)
      print(paste(now_row,"remain after collapsing on biopsy_id"))
    }

  }else{
    lost = data.frame()
  }
  expression_wider = filter(expression_wider,Available==TRUE) %>% 
    select(-Available) %>%
    mutate(sample_id = mrna_sample_id) %>% 
    unique()
  #lost = anti_join(original,expression_wider,
  #                 by=c("mrna_sample_id","seq_type","biopsy_id","protocol","ffpe_or_frozen")) %>% 
  #  unique()
  
  #join to other seq_types where possible
  capture_meta = get_gambl_metadata(seq_type_filter = "capture") %>% 
    filter(patient_id %in% expression_wider$patient_id) %>%
    select(sample_id,patient_id,biopsy_id) %>%
    rename(c("capture_sample_id"="sample_id"))
  expression_wider = left_join(expression_wider,capture_meta)
  
  genome_meta = get_gambl_metadata(seq_type_filter = "genome") %>% 
    filter(patient_id %in% expression_wider$patient_id) %>%
    select(sample_id,patient_id,biopsy_id) %>%
    rename(c("genome_sample_id"="sample_id"))
  expression_wider = left_join(expression_wider,genome_meta)
  
  return(select(expression_wider,-multi_exp))
}

#' @title Get Gene Expression.
#'
#' @description Get the expression for one or more genes for all GAMBL samples.
#'
#' @details Effectively get gene expression for one or multiple genes for all GAMBL samples.
#' For examples and more info, refer to the parameter descriptions as well as vignette examples.
#' 
#'
#' @param these_samples_metadata The data frame with sample metadata. Usually output of the get_gambl_metadata().
#' @param hugo_symbols One or more gene symbols.
#' @param ensembl_gene_ids One or more ensembl gene IDs. Only one of hugo_symbols or ensembl_gene_ids may be used.
#' @param all_genes Set to TRUE to return the full expression data frame without any subsetting. Avoid this if you don't want to use tons of RAM.
#'
#' @return A data frame with gene expression.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @rawNamespace import(vroom, except = c("col_skip", "fwf_positions", "default_locale", "date_names_lang", "cols_only", "output_column", "col_character", "col_guess", "spec", "as.col_spec", "fwf_cols", "cols", "col_date", "col_datetime", "locale", "col_time", "cols_condense", "col_logical", "col_number", "col_integer", "col_factor", "fwf_widths", "date_names_langs", "problems", "date_names", "col_double", "fwf_empty"))
#' @import dplyr readr tidyr GAMBLR.data
#' @export
#'
#' @examples
#' MYC_expr = get_gene_expression(hugo_symbols = "MYC", join_with = "mrna")
#'
#' #Read full expression values df (no subsetting on genes)
#' full_expression_df = get_gene_expression(all_genes = TRUE,
#'                                          join_with = "genome")
#'
#' #Use loaded df (in the previous step) to get expression values for IRF4 and MYC.
#' irf4_myc_expressions = get_gene_expression(hugo_symbols = c("IRF4", "MYC"),
#'                                            all_genes = FALSE,
#'                                            join_with = "genome",
#'                                            from_flatfile = FALSE,
#'                                            expression_data = full_expression_df)
#'
get_gene_expression = function(these_samples_metadata,
                               hugo_symbols,
                               ensembl_gene_ids,
                               all_genes = FALSE,
                               verbose=FALSE){
  
  if(!missing(hugo_symbols)){
    hugo_symbols = as.character(hugo_symbols)
  }else if(!missing(ensembl_gene_ids)){
    ensembl_gene_ids = as.character(ensembl_gene_ids)
  }
  
  sample_details = check_gene_expression() 
  # this contains all available non-redundant RNA-seq sample_ids. 
  # if necessary, subset it to samples in these_samples_metadata
  # The subsetting must consider the seq_type of each row in the metadata because there are as many as 3 sample IDs
  
  sample_details = filter(sample_details,
                          mrna_sample_id %in% these_samples_metadata$sample_id | 
                          capture_sample_id %in% these_samples_metadata$sample_id |
                          genome_sample_id %in% these_samples_metadata$sample_id)
  if(verbose){
    remaining_rows = nrow(sample_details)
    message(paste(remaining_rows,"samples from your metadata have RNA-seq data available"))
  }
  load_expression_by_samples = function(genes,samples,verbose){
      split_dir = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/DESeq2-0.0_salmon-1.0/mrna--gambl-icgc-all/split"
      all_sample_exp = list()
      for(sample in samples){
        sfile = paste0(split_dir,"/",sample,".tsv")
        if(!file.exists(sfile)){
          stop(paste("Can't find file for",sample))
        }
        rows = suppressMessages(read_tsv(sfile,progress = F,
                                         col_types = "cccncccc",
                                         num_threads = 8,
                                         lazy=TRUE)) %>% 
          filter(Hugo_Symbol %in% genes)
        all_sample_exp[[sample]] = rows
        nr = nrow(all_sample_exp[[sample]] )
        if(verbose){
          print(paste("loading data from",sample))
        }
      }
    all_rows = do.call("bind_rows",all_sample_exp)
    
    return(all_rows)
  }

  expression_long = load_expression_by_samples(genes=hugo_symbols,
                                               samples=sample_details$mrna_sample_id,
                                               verbose=verbose)
  
  expression_wide = pivot_wider(expression_long,
                                -ensembl_gene_id,
                                names_from="Hugo_Symbol",
                                values_from="expression")
  expression_wide = left_join(sample_details,expression_wide)
  
  return(expression_wide)
}