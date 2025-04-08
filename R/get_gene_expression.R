
#' @title Get Gene Expression.
#'
#' @description Get the expression for one or more genes for all GAMBL samples.
#'
#' @details Efficiently retrieve variance-stabilized and batch effect
#' corrected gene expression values for one, multiple or all genes
#' for all GAMBL samples.
#' For more information, refer to the parameter descriptions and examples.
#'
#' Warnings:
#'
#' 1) The speed of loading data is heavily impacted by how many samples
#' you load. For the sake of efficiency, be sure not to specify
#' extraneous samples.
#' 2) To reduce impact on memory (RAM), load only the data for the genes
#' you need.
#' 3) Combining lazy_join with all_genes will result in a data table
#' with samples on rows and genes on columns. *Use with caution*.
#' This is practically guaranteed to use more RAM than you want.
#' 4) Before you run this function, it's recommended that you run
#' `check_gene_expression` to determine which samples are available  
#'
#' @param these_samples_metadata The data frame with sample metadata.
#' Usually output of the get_gambl_metadata().
#' @param hugo_symbols One or more gene symbols.
#' Cannot be used in conjunction with ensembl_gene_ids. 
#' @param ensembl_gene_ids One or more ensembl gene IDs.
#' Cannot be used in conjunction with hugo_symbols. 
#' @param all_genes Set to TRUE for the full expression data without
#' any subsetting (*see warnings below*).
#' @param engine Either readr or grep. The grep engine usually will increase
#' the speed of loading but doesn't work if you want all genes or a very
#' long list.
#' @param format Either `wide` or `long`. Wide format returns one column
#' of expression values per gene. Long format returns one column of expression
#' values with the gene stored in a separate column. 
#' @param lazy_join If TRUE, your data frame will also have capture_sample_id
#' and genome_sample_id columns provided. See `check_gene_expression` for more
#' information.
#' @param arbitrarily_pick A stop-gap for handling the rare scenario where
#' the same Hugo_Symbol has more than one ensembl_gene_id. Set to TRUE only
#' if you encounter an error that states "Values are not uniquely identified;
#' output will contain list-cols."
#' @param HGNC When you request the wide matrix and all genes, this forces
#' the columns to contain hgnc_id rather than ensembl_gene_id
#' @param verbose Set to TRUE for a more chatty output
#' @param ... Optional parameters to pass along to `get_gambl_metadata` (only used in conjunction with lazy_join)
#'
#' @return A data frame with the first 9 columns identical to the columns from check_gene_expression and the remaining columns containing the expression values for each gene requested. 
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' # Get the expression for a single gene for every sample with RNA-seq
#' data in GAMBL
#' # This uses the default (grep) engine, which may be intolerably
#' slow on some systems
#' SOX11_exp_all = get_gene_expression(hugo_symbols = "SOX11")
#'                                        
#' # Get the expression for a few genes for all available samples AND get
#' # all available linkages to genome/capture samples without dropping anything
#' my_fave_gene_exp_long = get_gene_expression(
#'                         hugo_symbols = c("MYC","BCL2","EZH2"),
#'                         lazy_join=TRUE,
#'                         format="long")
#' 
#' # Get the expression values for the Wright gene set from every
#' # sample in the DLBCL_DLC cohort
#' # This is one example where the grep engine is significantly slower
#' than using the readr engine. This example shows the more efficient approach:
#' my_genes = GAMBLR.data::wright_genes_with_weights$Hugo_Symbol
#' my_meta = get_gambl_metadata() %>%
#'   dplyr::filter(cohort=="DLBCL_DLC"),
#'     seq_type %in% c('genome','capture')
#' wright_expr_with_DNA = get_gene_expression(engine="readr",
#'                                            hugo_symbols = my_genes,
#'                                            these_samples_metadata = my_meta)
#'
#' #Load the full expression values for every FL sample and all genes
#' # When tested on a gphost, this took less than a minute to run
#' my_meta = get_gambl_metadata() %>%
#'   dplyr::filter(pathology=="FL",
#'     seq_type=="mrna")
#' FL_expression_df = get_gene_expression(these_samples_metadata = my_meta,
#'                                        all_genes = TRUE)
#'
#' # Get the full expression table for every sample available
#' # in GAMBL (in the wide format)
#' all_exp_wide = get_gene_expression(all_genes=T)
#'
#' # Get the full expression table for every sample available in GAMBL
#' # (in the wide format) AND lazy-join to the minimal metadata
#' # NOTE: This will transpose the matrix, which makes it significantly slower.
#' # Also, due to incomplete Hugo_Symbols, it has to use the ENSG
#' # identifiers for column names. 
#' all_exp_wide = get_gene_expression(all_genes=T,lazy_join=T)
#'
#'
#' # If you want hgnc_symbol instead of Ensembl_gene_id you need to force
#' # the function to arbitrarily drop duplicates. *Not ideal*
#' exp_hgnc = get_gene_expression(these_samples_metadata = get_gambl_metadata(),
#'                                    all_genes = T,
#'                                    lazy_join = T,
#'                                    HGNC=TRUE,
#'                                    arbitrarily_pick = T)
#'
#' }
get_gene_expression = function(these_samples_metadata,
                               hugo_symbols,
                               ensembl_gene_ids,
                               all_genes = FALSE,
                               verbose=FALSE,
                               engine="grep",
                               format="wide",
                               lazy_join = FALSE,
                               arbitrarily_pick = FALSE,
                               HGNC = FALSE,
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
    if(format == "long") {
      stop("long format is only available when you supply a set of genes with ensembl_gene_ids or hugo_symbols")
    }
  }
  if(lazy_join){
    sample_details = check_gene_expression(show_linkages = T, ...) 
  }else{
    sample_details = check_gene_expression(show_linkages = F) 
  }
  if(verbose){
    message(paste("will attach these columns:",paste(colnames(sample_details),collapse=",")))
  }
  
  # this contains all available non-redundant RNA-seq sample_ids. 
  # if necessary, subset it to samples in these_samples_metadata
  # The subsetting must consider the seq_type of each row in the metadata because there are as many as 3 sample IDs
  
  if(!missing(these_samples_metadata)){
    original_row_num = nrow(sample_details)
    sample_details = inner_join(select(these_samples_metadata,biopsy_id,patient_id) %>% unique(),sample_details) 
    new_row_num = nrow(sample_details)
    if(verbose){
      print(paste("originally",original_row_num, "rows.", "After subsetting with the provided metadata,",new_row_num,"rows remain"))
    }
  }
  
  if(verbose){
    remaining_rows = nrow(sample_details)
    message(paste(remaining_rows,"samples from your metadata have RNA-seq data available"))
  }
  load_expression = function(genes,
                             id_type,
                             verbose,
                             engine=engine,
                             all=all_genes){
    tidy_expression_path = check_config_and_value("results_merged$tidy_expression_path")
    tidy_expression_path = str_remove(tidy_expression_path,".gz$")
    base_path = GAMBLR.helpers::check_config_and_value("project_base")
    #automatically default to file in tempfs if available
    tidy_expression_file = "/dev/shm/vst-matrix-Hugo_Symbol_tidy.tsv"
    if(file.exists(tidy_expression_file)){
      message(paste("using file in tempfs:",tidy_expression_file))
    }else{
      tidy_expression_file = paste0(base_path,tidy_expression_path)
    }
      if(engine=="grep"){
        if(id_type=="hugo" || id_type == "ensembl"){
          gene_ids = genes
        }else{
          stop("grep is only compatible with gene subsetting")
        }
        
       
        if(all_genes){
          stop("grep is only compatible with gene subsetting")
        }else{
          grep_cmd <- paste(gene_ids, collapse = " -e ") %>% 
            gettextf("grep -h -w -F -e Hugo_Symbol -e %s %s", . , tidy_expression_file)
          if(verbose){
            print(grep_cmd)
          }
          all_rows = fread(cmd = grep_cmd,verbose = F) 
          
        }
      }else{
        
        all_rows = read_tsv(tidy_expression_file,
                                         col_types = "cccncccc",
                                         num_threads = 8,
                                         lazy=TRUE)
      }
      if(id_type=="hugo"){
        if(verbose){
          message(paste("filtering to keep only these genes",paste(genes,collapse=",")))
        }
        all_rows = all_rows %>% 
            filter(Hugo_Symbol %in% genes)
      }else if(id_type=="ensembl"){
        if(verbose){
          message(paste("filtering to keep only these genes",paste(genes,collapse=",")))
        }
        all_rows = all_rows %>% 
            filter(ensembl_gene_id %in% genes)
      }
    return(all_rows)
  }
  if(!missing(hugo_symbols)){
    expression_long = load_expression(genes=hugo_symbols,
                                      id_type="hugo",
                                               verbose=verbose,
                                               engine=engine) %>%
      left_join(sample_details,.)
      if(arbitrarily_pick){
        if(verbose){
          message("if any hugo_symbols mapping to >1 ENSG are identified, the first one encountered will be arbitrarily used")
        }
        expression_wide = filter(expression_long,!is.na(expression)) %>% 
          group_by(Hugo_Symbol,mrna_sample_id,biopsy_id,patient_id) %>% 
          slice_head(n=1) %>%
          pivot_wider(.,
                      id_cols=-ensembl_gene_id,
                      names_from="Hugo_Symbol",
                      values_from="expression")
      }else{
        expression_wide = filter(expression_long,!is.na(expression)) %>% 
                               pivot_wider(.,
                                  id_cols=-ensembl_gene_id,
                                  names_from="Hugo_Symbol",
                                  values_from="expression")
      }
    
  }else if(!missing(ensembl_gene_ids)){
    expression_long = load_expression(genes=ensembl_gene_ids,
                                      id_type="ensembl",
                                                 verbose=verbose,
                                                 engine=engine) %>%
      left_join(sample_details,.)
    expression_wide = filter(expression_long,!is.na(expression)) %>% 
                               pivot_wider(.,
                                  id_cols=-Hugo_Symbol,
                                  names_from="ensembl_gene_id",
                                  values_from="expression")
  }else{ #all genes
    #just directly load the matrix, skipping the tidy format version and unnecessary pivoting
    base_path = GAMBLR.helpers::check_config_and_value("project_base")
    wide_expression_path = check_config_and_value("results_merged$ex_matrix_path")
    wide_expression_file = paste0(base_path,wide_expression_path)
    message(paste("loading all expression data from",wide_expression_file))
    expression_wide = suppressMessages(read_tsv(wide_expression_file))  %>%
      filter(!str_detect(gene_id, "PAR_Y"))
   
    #note that this takes up a lot less RAM than I was expecting:
    #format(object.size(all_exp),units="Gb")
    #[1] "1.2 Gb"
    nc = ncol(expression_wide) - 3
    message(paste("full matrix has",nc,"columns. Subsetting based on the prioritized results."))
    keep_ids = sample_details$mrna_sample_id
    expression_wide = expression_wide[,c(1,2,3,which(colnames(expression_wide) %in% keep_ids))]
    nc = ncol(expression_wide) - 3
    message(paste("kept",nc,"columns"))
    if(lazy_join){
      if(HGNC){
        message("transposing, setting hgnc_id as column name")
        if(arbitrarily_pick){
          expression_wide = expression_wide %>% 
            group_by(hgnc_symbol) %>% 
            slice_head(n=1) %>%
            ungroup() %>%
            filter(!is.na(hgnc_symbol))
        }
        expression_wide = select(expression_wide,-gene_id,-ensembl_gene_id) %>% 
          column_to_rownames("hgnc_symbol") %>%
          t() %>%
          as.data.table(keep.rownames=T) %>% 
          dplyr::rename("sample_id"="rn") %>%
          left_join(sample_details,.,by="sample_id")
      }else{
        message("transposing, setting ensembl_gene_id as column name")
        expression_wide = select(expression_wide,-gene_id,-hgnc_symbol) %>% 
          column_to_rownames("ensembl_gene_id") %>%
          t() %>%
          as.data.table(keep.rownames=T) %>% 
          dplyr::rename("sample_id"="rn") %>%
          left_join(sample_details,.,by="sample_id")
      }
      
      return(expression_wide)
    }else{
        return(expression_wide)
    }
  }
  if(format == "long") {
    return(expression_long)
  }else{
    return(expression_wide)
  }
}