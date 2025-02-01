
#' @title Get Raw Expression Counts
#'
#' @description Get the raw read counts from RNA-seq for one or more genes for all GAMBL samples.
#'
#' @details Efficiently retrieve raw gene expression values (read counts) for one, multiple or all genes for all GAMBL samples.
#' For examples and more info, refer to the parameter descriptions as well as vignette examples.
#'
#' @param these_samples_metadata The data frame with sample metadata. Usually output of the get_gambl_metadata().
#' @param sample_id_column Specify which column in your metadata contains the sample_id you want used instead of the mrna sample_id
#' @param from_flatfile Set to FALSE to use the database instead of reading from flatfiles
#' @param check For basic debugging. Set to TRUE to obtain basic information about the number of samples in your metadata with expression data available 
#' @param all_samples Set to TRUE to force the function to return all available data (should rarely be necessary)
#' @param map_to_symbol Set to TRUE to obtain the mappings between the rows in the count matrix and HGNC gene symbol/alias
#' @param TPM Set to TRUE to get TPM estimates instead of counts
#' 
#' @return A list containing a counts matrix and the associated metadata for DESeq2
#'
#' @import dplyr readr tidyr msigdbr org.Hs.eg.db
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
                                     existing_sample_id_column,
                                     new_sample_id_column,
                                     all_samples=FALSE,
                                     check=FALSE,
                                     from_flatfile=TRUE,
                                     map_to_symbol=FALSE,
                                     verbose=FALSE,
                                     TPM = FALSE){
  if(missing(these_samples_metadata)){
    
    if(!all_samples & !check){
      warning("Missing these_samples_metadata. Results will contain data from all available samples. ")
      warning("This is almost certainly NOT what you want!")
      stop("re-run this function with all_samples = TRUE if you wish to proceed")
    }
    these_samples_metadata = get_gambl_metadata() %>% 
      dplyr::filter(seq_type == "mrna")
    
  }
  if(!missing(existing_sample_id_column)){
    sample_ids = pull(these_samples_metadata,{{existing_sample_id_column}}) 
  }else{
    sample_ids = pull(these_samples_metadata,sample_id)
  }
  
  
  if(from_flatfile){
    files_dir = paste0(config::get("project_base"),
                       GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$expression$salmon$counts))
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
    if(verbose){
      #return(file_df)
      print(tail(file_df,150))
    }
    
    df_list <- suppressMessages(map(file_df$full_path, read_tsv))
    expression_long = bind_rows(df_list)
    if(TPM){
      expression_long = expression_long %>% dplyr::select(-count)
      message("transposing")
      expression_wide = pivot_wider(expression_long,
                                    names_from="sample_id",
                                    values_from="TPM") %>%
        column_to_rownames("gene")
    }else{
      expression_long = expression_long %>% dplyr::select(-TPM)
      
      if(!missing(new_sample_id_column)){

        #swap the sample_id before transposing
        if(any(duplicated(these_samples_metadata[,new_sample_id_column]))){
          stop("new sample ids must all be unique")
        }
        message("ID substitution")
        expression_long = left_join(expression_long,select(these_samples_metadata,
                                                           {{existing_sample_id_column}},
                                                           {{new_sample_id_column}}),
                                    by=c("sample_id"={{existing_sample_id_column}})) %>%
          select(-sample_id)
        message("transposing")
        expression_wide = pivot_wider(expression_long,
                                      names_from=new_sample_id_column,
                                      values_from="count") %>%
          column_to_rownames("gene")

        expression_metadata =  dplyr::filter(these_samples_metadata, 
                                             .data[[new_sample_id_column]] %in% colnames(expression_wide)) %>%
          mutate(sample_id=.data[[new_sample_id_column]] ) %>%
          column_to_rownames({{new_sample_id_column}})
        expression_metadata = expression_metadata[colnames(expression_wide),]
      }else{
        expression_wide = pivot_wider(expression_long,
                                      names_from="sample_id",
                                      values_from="count") %>%
          column_to_rownames("gene")
        expression_metadata = these_samples_metadata %>%
          dplyr::filter(sample_id %in% colnames(expression_wide)) %>%
          column_to_rownames("sample_id")
        #expression_metadata = expression_metadata[colnames(expression_wide),]
      }
      
      
    }
    
    
    
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
  
  if(map_to_symbol){
    ensembl_id_clean = gsub("\\..+","",rownames(expression_wide))
    gene_symbol = mapIds(org.Hs.eg.db,
                         keys=ensembl_id_clean,
                         keytype="ENSEMBL",
                         column="SYMBOL")
    mapped = data.frame(original_id=rownames(expression_wide),
                        clean_id=names(gene_symbol),
                        symbol=unname(gene_symbol))
    return(list(counts=expression_wide,metadata=expression_metadata,IDs=mapped))
  }else{
    return(list(counts=expression_wide,metadata=expression_metadata))
  }
  
  
}
