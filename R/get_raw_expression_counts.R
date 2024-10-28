
#' @title Get Raw Expression Counts
#'
#' @description Get the raw read counts from RNA-seq for one or more genes for all GAMBL samples.
#'
#' @details Efficiently retrieve raw gene expression values (read counts) for one, multiple or all genes for all GAMBL samples.
#' For examples and more info, refer to the parameter descriptions as well as vignette examples.
#'
#' @param these_samples_metadata The data frame with sample metadata. Usually output of the get_gambl_metadata().
#' @param from_flatfile Set to FALSE to use the database instead of reading from flatfiles
#' @param check For basic debugging. Set to TRUE to obtain basic information about the number of samples in your metadata with expression data available 
#' @param all_samples Set to TRUE to force the function to return all available data (should rarely be necessary)
#' 
#' @return A list containing a counts matrix and the associated metadata for DESeq2
#'
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' 
#' schmitz_meta = get_gambl_metadata() %>% 
#'     filter(seq_type=="mrna",cohort=="dlbcl_schmitz")
#' exp_out = get_raw_expression_counts(these_samples_metadata = schmitz_meta)
#' 
#' # Create DESeq data set directly from the two named objects in the output
#' 
#' dds <- DESeqDataSetFromMatrix(countData = exp_out$counts,
#'     colData = exp_out$metadata,
#'     design = ~ COO_consensus + sex)
#'     
#' # Run a basic DESeq analysis
#' dds <- DESeq(dds)
#' res <- results(dds, 
#'     name="COO_consensus_GCB_vs_ABC",
#'     lfcThreshold=2,alpha=0.1)
#' # Filter outputs using padj, logFC and baseMean (more highly expressed overall)     
#' res_df = as.data.frame(res) %>% 
#'     filter(padj<0.1,baseMean>500)
#'
#' show_genes = rownames(res_df)
#' vsd <- vst(dds, blind=FALSE)
#' 
#' #Visualize the results with a heatmap
#' column_ha = HeatmapAnnotation(df=select(exp_out$metadata,COO_consensus,sex))
#' Heatmap(assay(vsd)[show_genes,],
#'     row_names_gp = gpar(fontsize=5),
#'     bottom_annotation = column_ha,
#'     show_column_names = F)
#' 
#'
get_raw_expression_counts = function(these_samples_metadata,
                                     all_samples=FALSE,
                                     check=FALSE,
                                     from_flatfile=FALSE){
  if(missing(these_samples_metadata)){
    
    if(!all_samples & !check){
      warning("Missing these_samples_metadata. Results will contain data from all available samples. ")
      warning("This is almost certainly NOT what you want!")
      stop("re-run this function with all_samples = TRUE if you wish to proceed")
    }
    these_samples_metadata = get_gambl_metadata() %>% 
      dplyr::filter(seq_type == "mrna")
    
  }
  sample_ids = pull(these_samples_metadata,sample_id)
  
  if(from_flatfile){
    files_dir = paste0(config::get("project_base"),GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$expression$salmon$counts))
    all_files = dir(files_dir)
    file_df = data.frame(name=all_files) %>% 
      mutate(sample_id = str_remove(name,".tsv.gz"))
    if(check){
      all_n = nrow(file_df)
      file_df = filter(file_df,sample_id %in% sample_ids) 
      these_n = nrow(file_df)
      print(paste("TOTAL AVAILABLE:",all_n))
      print(paste("samples matching provided metadata:",these_n))
      return()
    }
    file_df = filter(file_df,sample_id %in% sample_ids) 
    #load and combine
    file_df = mutate(file_df,full_path = paste0(files_dir,name))
    df_list <- map(file_df$full_path, read_tsv(show_col_types=F))
    expression_long = bind_rows(df_list)
    message("transposing")
    expression_wide = pivot_wider(expression_long,names_from="sample_id",values_from="count") %>%
      column_to_rownames("gene")
    expression_metadata = dplyr::filter(these_samples_metadata,sample_id %in% colnames(expression_wide)) %>%
      column_to_rownames("sample_id")
    expression_metadata = expression_metadata[colnames(expression_wide),]
    
  }else{
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = "gambl_test")
    counts_table <- tbl(con, "salmon_counts")
    if(check){
      expression_one = dplyr::filter(counts_table,gene=="ENSG00000001084.13") %>% as.data.frame()
      all_n = nrow(expression_one)
      expression_one = filter(expression_one,sample_id %in% sample_ids) 
      these_n = nrow(expression_one)
      print(paste("TOTAL AVAILABLE:",all_n,"MATCHING METADATA:",these_n))
    }else{
      message("querying database")
      sample_ids = sample_ids[c(1:1000)]
      expression_long = dplyr::filter(counts_table,
                                      sample_id %in% sample_ids) %>%
        as.data.frame()
      message("transposing")
      expression_wide = pivot_wider(expression_long,names_from="sample_id",values_from="count") %>%
        column_to_rownames("gene")
      expression_metadata = dplyr::filter(these_samples_metadata,sample_id %in% colnames(expression_wide)) %>%
        column_to_rownames("sample_id")
      expression_metadata = expression_metadata[colnames(expression_wide),]
    }
    
    
    dbDisconnect(con)
  }
  
  
  
  return(list(counts=expression_wide,metadata=expression_metadata))
}
