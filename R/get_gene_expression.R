
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
#' @import dplyr readr tidyr
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
                               engine="grep",
                               from_merge=TRUE){
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
      #split_dir = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/DESeq2-0.0_salmon-1.0/mrna--gambl-icgc-all/split/"
      if(engine=="grep"){
        if(!missing(hugo_symbols)){
          gene_ids = hugo_symbols
        }else if(!missing(ensembl_gene_ids)){
          gene_ids = ensembl_gene_ids
        }else{
          stop("grep is only compatible with gene subsetting")
        }
        sample_files = paste(paste0(split_dir,samples,".tsv"),collapse=" " )
        sample_file = paste(paste0(split_dir,samples[1],".tsv"),collapse=" " )
        
        if(from_merge){
          
          tidy_expression_path = check_config_value(config::get("results_merged")$tidy_expression_path)
          tidy_expression_path = str_remove(tidy_expression_path,".gz$")
          base_path = GAMBLR.helpers::check_config_value(config::get("project_base"))
          tidy_expression_file = paste0(base_path,tidy_expression_path)
          grep_cmd <- paste(gene_ids, collapse = " -e ") %>% 
            gettextf("grep -h -w -F -e Hugo_Symbol -e %s %s", . , tidy_expression_file)
          if(verbose){
            print(grep_cmd)
          }
          all_rows = fread(cmd = grep_cmd,verbose = T) 
        }else{
          all_sample_exp = list()
          for(sample in samples){
            sfile = paste0(split_dir,"/",sample,".tsv")
            grep_cmd <- paste(gene_ids, collapse = " -e ") %>% 
              gettextf("grep -h -w -e %s %s", . , sfile)
            if(verbose){
              #print(grep_cmd)
              print(paste("loading data from",sample))
            }
            rows = fread(cmd=grep_cmd,col.names = c("ensembl_gene_id",
                                                    "Hugo_Symbol",     
                                                    "mrna_sample_id",  
                                                    "expression",      
                                                    "patient_id",     
                                                    "biopsy_id",
                                                    "protocol",
                                                    "ffpe_or_frozen"))
            
            all_sample_exp[[sample]] = rows
            nr = nrow(all_sample_exp[[sample]] )
            
          }
          all_rows = do.call("bind_rows",all_sample_exp)
          #grep_cmd <- paste(gene_ids, collapse = " -e ") %>% 
          #  gettextf("grep -h -w -F -e %s %s)",sample_file,., sample_files)
          #if(verbose){
          #  print(grep_cmd)
          #}
          #all_rows = fread(cmd = grep_cmd,verbose=T)
        }
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