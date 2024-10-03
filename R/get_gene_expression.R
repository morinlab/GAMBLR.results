
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
#' 3) Before you run this function, it's recommended that you run `check_gene_expression` to determine which samples are available  
#'
#' @param these_samples_metadata The data frame with sample metadata. Usually output of the get_gambl_metadata().
#' @param hugo_symbols One or more gene symbols. Cannot be used in conjunction with ensembl_gene_ids. 
#' @param ensembl_gene_ids One or more ensembl gene IDs. Cannot be used in conjunction with hugo_symbols. 
#' @param all_genes Set to TRUE to return the full expression data frame without any subsetting (see warnings below). 
#' @param engine Either readr or grep. The grep engine usually will increase the speed of loading but doesn't work if you want all genes or a very long list.
#' @param lazy_join If TRUE, your data frame will also have capture_sample_id and genome_sample_id columns provided. See `check_gene_expression` for more information.
#' @param ... Optional parameters to pass along to `get_gambl_metadata` (only used in conjunction with lazy_join)
#'
#' @return A data frame with the first 9 columns identical to the columns from check_gene_expression and the remaining columns containing the expression values for each gene requested. 
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' 
#' # Get the expression for a single gene for every sample with RNA-seq data in GAMBL
#' # When tested on a gphost, this took about 45 minutes to run with the grep engine
#' SOX11_exp_all = get_gene_expression(hugo_symbols = "SOX11")
#'                                        
#' Get the expression for a single gene for every sample wit RNA-seq data AND get all available linkages to genome/capture samples without dropping anything
#' SOX11_exp_all = get_gene_expression(hugo_symbols = "SOX11",lazy_join=TRUE)
#' 
#' # Get the expression values for the Wright gene set from every DLBCL sample with either genome or capture data in GAMBL.
#' wright_gene_expr_all_DLBCL_with_DNA = get_gene_expression(hugo_symbols = GAMBLR.data::wright_genes_with_weights$Hugo_Symbol,
#'                                                  these_samples_metadata = get_gambl_metadata() %>% 
#'                                                  dplyr::filter(pathology=="DLBCL"),seq_type %in% c('genome','capture'))
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
                               lazy_join = FALSE,
                               ...){
  if(missing(these_samples_metadata)){
    warning("Missing these_samples_metadata. Results will contain data from all available samples. ")
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
  if(lazy_join){
    sample_details = check_gene_expression(show_linkages = T, ...) 
  }else{
    sample_details = check_gene_expression(show_linkages = F) 
  }
  
  # this contains all available non-redundant RNA-seq sample_ids. 
  # if necessary, subset it to samples in these_samples_metadata
  # The subsetting must consider the seq_type of each row in the metadata because there are as many as 3 sample IDs
  
  if(!missing(these_samples_metadata)){
    original_row_num = nrow(sample_details)
    sample_details = inner_join(select(these_samples_metadata,biopsy_id,patient_id) %>% unique(),sample_details) 
    #filter(sample_details,
    #                        mrna_sample_id %in% these_samples_metadata$sample_id | 
    #                          capture_sample_id %in% these_samples_metadata$sample_id |
    #                          genome_sample_id %in% these_samples_metadata$sample_id)
    new_row_num = nrow(sample_details)
    if(verbose){
      print(paste("originally",original_row_num, "rows.", "After subsetting with the provided metadata,",new_row_num,"rows remain"))
    }
  }
  
  if(verbose){
    remaining_rows = nrow(sample_details)
    message(paste(remaining_rows,"samples from your metadata have RNA-seq data available"))
  }
  load_expression_by_samples = function(samples,hugo_symbols,ensembl_gene_ids,verbose,read_engine=engine){
    tidy_expression_path = check_config_value(config::get("results_merged")$tidy_expression_path)
    base_path = GAMBLR.helpers::check_config_value(config::get("project_base"))
    tidy_expression_file = paste0(base_path,tidy_expression_path)
    print(tidy_expression_file)
      if(read_engine=="grep"){
        if(!missing(hugo_symbols)){
          gene_ids = hugo_symbols
        }else if(!missing(ensembl_gene_ids)){
          gene_ids = ensembl_gene_ids
        }else{
          stop("grep is only compatible with gene subsetting")
        }
        
        
        grep_cmd <- paste(gene_ids, collapse = " -e ") %>% 
          gettextf("grep -h -w -F -e Hugo_Symbol -e %s %s", . , tidy_expression_file)
        if(verbose){
          print(grep_cmd)
        }
        all_rows = fread(cmd = grep_cmd,verbose = F) 
        
      }else{
        
        all_rows = read_tsv(tidy_expression_file,
                                         col_types = "cccncccc",
                                         num_threads = 8,
                                         lazy=TRUE)
          if(!missing(hugo_symbols)){
            all_rows = all_rows %>% 
              filter(Hugo_Symbol %in% hugo_symbols)
          }else if(!missing(ensembl_gene_ids)){
            all_rows = all_rows %>% 
              filter(ensembl_gene_id %in% ensembl_gene_ids)
          }

      }
    
    return(all_rows)
  }
  if(!missing(hugo_symbols)){
    expression_long = load_expression_by_samples(sample_details$mrna_sample_id,
                                                 hugo_symbols=hugo_symbols,
                                               verbose=verbose,
                                               read_engine=engine)
    check_config_value(config::get("results_merged")$tidy_expression_path)
    expression_wide = pivot_wider(expression_long,
                                  -ensembl_gene_id,
                                  names_from="Hugo_Symbol",
                                  values_from="expression")
  }else if(!missing(ensembl_gene_ids)){
    expression_long = load_expression_by_samples(sample_details$mrna_sample_id,
                                                 ensembl_gene_ids=ensembl_gene_ids,
                                
                                                 verbose=verbose,
                                                 read_engine=engine)
    expression_wide = pivot_wider(expression_long,
                                  -Hugo_Symbol,
                                  names_from="ensembl_gene_id",
                                  values_from="expression")
  }else{
    expression_long = load_expression_by_samples(sample_details$mrna_sample_id)
    #expression_wide = pivot_wider(expression_long,
    #                              -Hugo_Symbol,
    #                              names_from="ensembl_gene_id",
    #                              values_from="expression")
    
    expression_wide = pivot_wider(expression_long,
                                  -ensembl_gene_id,
                                  names_from="Hugo_Symbol",
                                  values_from="expression")
  }
  
  
  expression_wide = left_join(sample_details,expression_wide)
  
  return(expression_wide)
}