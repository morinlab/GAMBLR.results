
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

#' @title Get Gene Expression.
#'
#' @description Get the expression for one or more genes for all GAMBL samples.
#'
#' @details Efficiently retrieve variance-stabilized and batch effect corrected gene expression values for one, multiple or all genes for all GAMBL samples.
#' For examples and more info, refer to the parameter descriptions as well as vignette examples.
#'
#' Warnings: 
#'
#' 1) The speed of loading data is heavily impacted by how many samples you load. For the sake of efficiency, be sure not to specify extraneous samples. 
#' 2) To reduce impact on memory (RAM), load only the data for the genes you need. 
#' 3) To get the most complete set of RNA-seq data, be sure to call get_gambl_metadata(seq_type = "mrna") otherwise it will default to genomes. 
#' 4) Before you run this function, it's recommended that you run `check_gene_expression` to determine which samples are available  
#'
#' @param these_samples_metadata The data frame with sample metadata. Usually output of the get_gambl_metadata().
#' @param hugo_symbols One or more gene symbols. Cannot be used in conjunction with ensembl_gene_ids. 
#' @param ensembl_gene_ids One or more ensembl gene IDs. Cannot be used in conjunction with hugo_symbols. 
#' @param all_genes Set to TRUE to return the full expression data frame without any subsetting (see warnings below). 
#' @param engine Either readr or grep. The grep engine usually will increase the speed of loading but doesn't work if you want all genes or a very long list.
#'
#' @return A data frame with the first 9 columns identical to the columns from check_gene_expression and the remaining columns containing the expression values for each gene requested. 
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @rawNamespace import(vroom, except = c("col_skip", "fwf_positions", "default_locale", "date_names_lang", "cols_only", "output_column", "col_character", "col_guess", "spec", "as.col_spec", "fwf_cols", "cols", "col_date", "col_datetime", "locale", "col_time", "cols_condense", "col_logical", "col_number", "col_integer", "col_factor", "fwf_widths", "date_names_langs", "problems", "date_names", "col_double", "fwf_empty"))
#' @import dplyr readr tidyr GAMBLR.data
#' @export
#'
#' @examples
#' 
#' # Get the expression for a single gene for every sample with RNA-seq data in GAMBL
#' # When tested on a gphost, this took about 4 minutes to run with the readr engine
#' SOX11_exp_all = get_gene_expression(these_samples_metadata = 
#'                                        get_gambl_metadata(seq_type_filter="mrna"),
#'                                        hugo_symbols = "SOX11")
#'                                        
#' # Get the expression for a single gene for every sample with both RNA-seq data AND genome data in GAMBL
#' # When tested on a gphost, this took about 2 minutes to run with the default readr engine and 11 seconds with the grep engine
#' SOX11_exp_all = get_gene_expression(these_samples_metadata = 
#'                                        get_gambl_metadata(seq_type_filter="genome"),
#'                                        hugo_symbols = "SOX11",
#'                                        engine = "grep")
#' 
#' # Get the expression values for the Wright gene set from every DLBCL sample with either genome or capture data in GAMBL.
#' # When tested on a gphost, this took about 3 minutes to run with the readr engine and about a minute with the grep engine
#' wright_gene_expr_all_DLBCL_with_DNA = get_gene_expression(hugo_symbols = GAMBLR.data::wright_genes_with_weights$Hugo_Symbol,
#'                                                  these_samples_metadata = get_gambl_metadata(seq_type_filter = c('genome','capture')) %>% 
#'                                                  dplyr::filter(pathology=="DLBCL"))
#'
#' #Load the full expression values for every FL sample with no subsetting on genes
#' # When tested on a gphost, this took about 6 minutes to run
#' FL_expression_df = get_gene_expression(these_samples_metadata = 
#'                                                  get_gambl_metadata(seq_type_filter="mrna") %>% 
#'                                                      dplyr::filter(pathology=="FL"),
#'                                                  all_genes = TRUE)
#'
#'
get_gene_expression = function(these_samples_metadata,
                               hugo_symbols,
                               ensembl_gene_ids,
                               all_genes = FALSE,
                               verbose=FALSE,
                               engine="readr"){
  if(missing(these_samples_metadata)){
    stop("Missing these_samples_metadata. You must supply a data frame of metadata to specify which samples you want to work with. ")
  }
  if(!missing(hugo_symbols)){
    hugo_symbols = as.character(hugo_symbols)
  }else if(!missing(ensembl_gene_ids)){
    ensembl_gene_ids = as.character(ensembl_gene_ids)
  }else{
    #both are missing
    if(!all_genes){
      stop("You must either provide a vector of ensembl_gene_ids, hugo_symbols or explicitly set all_genes = TRUE (if you really want to load everything)")
    }
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
  load_expression_by_samples = function(hugo_symbols,ensembl_gene_ids,samples,verbose,engine=engine){
      split_dir = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/DESeq2-0.0_salmon-1.0/mrna--gambl-icgc-all/split/"
      if(engine=="grep"){
        if(!missing(hugo_symbols)){
          gene_ids = hugo_symbols
        }else if(!missing(ensembl_gene_ids)){
          gene_ids = ensembl_gene_ids
        }else{
          stop("grep is only compatible with gene subsetting")
        }
        sample_files = paste(paste0(split_dir,samples,".tsv"),collapse=" " )
        grep_cmd <- paste(gene_ids, collapse = " -e ") %>% 
          gettextf("grep -h -w -F -e %s %s", ., sample_files)
        if(verbose){
          print(grep_cmd)
        }
        
        all_rows = fread(cmd = grep_cmd,col.names = c("ensembl_gene_id",
                                                                  "Hugo_Symbol",
                                                                  "mrna_sample_id",
                                                                  "expression",
                                                                  "patient_id",
                                                                  "biopsy_id",
                                                                  "protocol",
                                                                  "ffpe_or_frozen"))
      }else{
        all_sample_exp = list()
        for(sample in samples){
          sfile = paste0(split_dir,"/",sample,".tsv")
          if(!file.exists(sfile)){
            stop(paste("Can't find file for",sample))
          }
          rows = suppressMessages(read_tsv(sfile,progress = F,
                                         col_types = "cccncccc",
                                         num_threads = 8,
                                         lazy=TRUE))
          if(!missing(hugo_symbols)){
            rows = rows %>% 
              filter(Hugo_Symbol %in% hugo_symbols)
          }else if(!missing(ensembl_gene_ids)){
            rows = rows %>% 
              filter(ensembl_gene_id %in% ensembl_gene_ids)
          }
          all_sample_exp[[sample]] = rows
          nr = nrow(all_sample_exp[[sample]] )
          if(verbose){
            print(paste("loading data from",sample))
          }
        }
        all_rows = do.call("bind_rows",all_sample_exp)
      }
    
    return(all_rows)
  }
  if(!missing(hugo_symbols)){
    expression_long = load_expression_by_samples(hugo_symbols=hugo_symbols,
                                               samples=sample_details$mrna_sample_id,
                                               verbose=verbose,
                                               engine=engine)
    expression_wide = pivot_wider(expression_long,
                                  -ensembl_gene_id,
                                  names_from="Hugo_Symbol",
                                  values_from="expression")
  }else if(!missing(ensembl_gene_ids)){
    expression_long = load_expression_by_samples(ensembl_gene_ids=ensembl_gene_ids,
                                                 samples=sample_details$mrna_sample_id,
                                                 verbose=verbose,
                                                 engine=engine)
    expression_wide = pivot_wider(expression_long,
                                  -Hugo_Symbol,
                                  names_from="ensembl_gene_id",
                                  values_from="expression")
  }else{
    expression_long = load_expression_by_samples(samples=sample_details$mrna_sample_id,
                                                 verbose=verbose)
    expression_wide = pivot_wider(expression_long,
                                  -Hugo_Symbol,
                                  names_from="ensembl_gene_id",
                                  values_from="expression")
  }
  
  
  expression_wide = left_join(sample_details,expression_wide)
  
  return(expression_wide)
}