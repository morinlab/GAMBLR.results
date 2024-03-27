#' @title Get Gene Expression.
#'
#' @description Get the expression for one or more genes for all GAMBL samples.
#'
#' @details Effectively get gene expression for one or multiple genes for all GAMBL samples.
#' This function can also take an already loaded expression matrix (`expression_data`)
#' to prevent the user from having to load the full expression matrix if this function needs to be run in an interactive session.
#' For examples and more info, refer to the parameter descriptions as well as vignette examples.
#' The function has argument `engine`, which accepts string "read_tsv", "grep", "vroom", and "fread". This will determine the way
#' the data is imported into R. When testing on GSC, the grep was the fastest but with a lot of variation in the run time (anywhere between 4-10 min).
#' Other engines produced similar run times (~ 7 min) on GSC, with vroom engine being the most consistent one. However, on other
#' systems (especially with fast hard drives and MacBooks for remote users) the read_tsv engine was significantly faster than the grep.
#' 
#' If `these_metadata_samples` is not provided, and `join_with` is one of `"mrna"`, `"genome"` or `"capture"`, 
#' `get_gambl_metadata` is called internally to retrieve all samples IDs of the respective seq type. If similar 
#' scenario but `join_with = NULL`, `get_gambl_metadata` is called to retrieve all sample IDs of all seq types.
#' For each metadata sample, `get_gene_expression` tries to return the expression of the genes provided by the 
#' `hugo_symbols` parameter. If a gene expression can not be retrieved for a sample, the function returns NA.
#' 
#' If `join_with = "mrna"`, `get_gene_expression` retrieves the gene expression of sample IDs from the metadata 
#' by directly matching them to mRNA sample IDs from the internal gene expression file. If `join_with` is other 
#' than "mrna" (genome", "capture", or NULL), the function links the sample IDs by matching both patient and 
#' biopsy IDs from the metadata and from the expression files.
#' 
#' The `prioritize_rows_by` parameter may be used to prioritize rows and avoid duplications in the output table. 
#' A duplication is when a same sample ID from the metadata is linked to more than one mRNA sample ID from the 
#' internal gene expression file, hence the metadata sample ID is associated to more than one different expression 
#' level. To filter out duplications, provide to `prioritize_rows_by` a named list of vectors, where a name 
#' specifies a column (contained in the output) and its respective vector elements refer to possible values of 
#' this column to be prioritized. The first values of the vector have higher prioritization. First, filtering is 
#' applied using the column specified by the first element of list `prioritize_rows_by`. If any duplication remains, 
#' the next element is used, and so on. This parameter is optional and if not provided, the filtering is not applied. 
#' If a duplication can not be solved, their rows will be marked as `1` in the output `multi_exp` column. Below is 
#' an example of how to set up the `prioritize_rows_by` parameter:
#' 
#' ```
#' prioritize_rows_by = list(
#'   protocol = c("Strand_Specific_Transcriptome_2",
#'                "Strand_Specific_Transcriptome_3"),
#'   ffpe_or_frozen = "frozen"
#' )
#' ```
#'
#' @param these_samples_metadata The data frame with sample metadata. Usually output of the get_gambl_metadata().
#' @param hugo_symbols One or more gene symbols.
#' @param ensembl_gene_ids One or more ensembl gene IDs. Only one of hugo_symbols or ensembl_gene_ids may be used.
#' @param engine Specific way to import the data into R. Defaults to "read_tsv". Other acceptable options are "grep", "vroom", and "fread".
#' @param join_with The seq type used to join the expression data to the metadata table. Can be one of NULL (default), 
#'   "mrna", "genome", or "capture". See the **Details** section for more information. 
#' @param all_genes Set to TRUE to return the full expression data frame without any subsetting. Avoid this if you don't want to use tons of RAM.
#' @param expression_data Optional argument to use an already loaded expression data frame (prevent function to re-load full df from flat file or database).
#' @param from_flatfile Deprecated but left here for backwards compatibility.
#' @param prioritize_rows_by A named list with one or more vectors. Provide this parameter if you want to filter out 
#'   duplications (as indicated by the `multi_exp` column of the output table) by prioritizing rows based on the values 
#'   of columns specified by this parameter. See the **Details** section for more information. 
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
                               engine = "read_tsv",
                               join_with = NULL,
                               all_genes = FALSE,
                               expression_data,
                               from_flatfile = TRUE,
                               prioritize_rows_by){
  
  # check parameters
  if(!is.null(join_with)){
    stopifnot("`join_with` must be one of NULL, \"mrna\", \"genome\", or \"capture\"." = 
                join_with %in% c("mrna", "genome", "capture"))
  }
  
  if(!missing(hugo_symbols)){
    hugo_symbols = as.character(hugo_symbols)
  }else if(!missing(ensembl_gene_ids)){
    ensembl_gene_ids = as.character(ensembl_gene_ids)
  }
  
  database_name = GAMBLR.helpers::check_config_value(config::get("database_name"))
  if(missing(these_samples_metadata)){
    if(is.null(join_with)){
      these_samples_metadata = get_gambl_metadata(seq_type_filter = c("genome", "capture", "mrna"), 
                                                  only_available = FALSE)
      
    }else if(join_with == "mrna"){
      these_samples_metadata = get_gambl_metadata(seq_type_filter = "mrna", only_available = FALSE)
      
    }else if(join_with == "genome"){
      these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome", only_available = FALSE)
      
    }else if(join_with == "capture"){
      these_samples_metadata = get_gambl_metadata(seq_type_filter = "capture", only_available = FALSE)
    }
  }
  
  if(missing(hugo_symbols) & missing(ensembl_gene_ids) & !all_genes){
    stop("ERROR: supply at least one gene symbol or Ensembl gene ID")
  }else if(!missing(hugo_symbols) & !missing(ensembl_gene_ids)){
    stop("ERROR: Both hugo_symbols and ensembl_gene_ids were provided. Please provide only one type of ID.")
  }
  #tidy_expression_file = config::get("results_merged")$tidy_expression_file
  #use combination of base path and relative path instead of full path for flexibility accross sites
  tidy_expression_path = GAMBLR.helpers::check_config_value(config::get("results_merged")$tidy_expression_path)
  base_path = GAMBLR.helpers::check_config_value(config::get("project_base"))
  tidy_expression_file = paste0(base_path,tidy_expression_path)
  tidy_expression_file = gsub(".gz$","",tidy_expression_file)
  
  #check permission and updates paths accordingly
  permissions = file.access(tidy_expression_file, 4)
  if(permissions == -1 ){
    message("restricting to non-ICGC data")
    tidy_expression_path = GAMBLR.helpers::check_config_value(config::get("results_merged")$tidy_expression_path_gambl)
    tidy_expression_file = paste0(base_path, tidy_expression_path)
  }
  
  if(!missing(expression_data)){
    tidy_expression_data = as.data.frame(expression_data) #is this necessary? Will it unnecessarily duplicate a large object if it's already a data frame?
    if(!missing(hugo_symbols)){
      #lazily filter on the fly to conserve RAM
      wide_expression_data = tidy_expression_data %>%
        dplyr::filter(Hugo_Symbol %in% hugo_symbols) %>%
        dplyr::select(-ensembl_gene_id) %>%
        group_by(mrna_sample_id,Hugo_Symbol) %>% #deal with non 1:1 mapping of Hugo to Ensembl
        slice_head() %>%
        as.data.frame() %>%
        pivot_wider(names_from = Hugo_Symbol, values_from = expression)
    }else if(!missing(ensembl_gene_ids)){
      wide_expression_data = tidy_expression_data %>%
        dplyr::filter(ensembl_gene_id %in% ensembl_gene_ids) %>%
        dplyr::select(-Hugo_Symbol) %>%
        as.data.frame() %>%
        pivot_wider(names_from = ensembl_gene_id, values_from = expression)
    }else{
      
      #for when a user wants everything. Need to handle the option of getting back Hugo_Symbol instead
      wide_expression_data = tidy_expression_data %>%
        dplyr::select(-Hugo_Symbol) %>%
        as.data.frame() %>%
        pivot_wider(names_from = ensembl_gene_id, values_from = expression)
    }
  }else{
    if(!file.exists(tidy_expression_file)){
      print(paste("missing: ", tidy_expression_file))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE= "remote")')
      check_host()
    }
    #only ever load the full data frame when absolutely necessary
    if(all_genes & missing(ensembl_gene_ids) & missing(hugo_symbols)){
      wide_expression_data = suppressMessages(read_tsv(tidy_expression_file)) %>%
        as.data.frame() %>%
        pivot_wider(names_from = ensembl_gene_id, values_from = expression)
    }else{
      if(!missing(hugo_symbols)){
        if(engine == "read_tsv"){
          message("Will read the data using read_tsv")
          wide_expression_data = read_tsv(tidy_expression_file,lazy=TRUE) %>%
            dplyr::filter(Hugo_Symbol %in% hugo_symbols)
        } else if(engine == "vroom"){
          message("Will read the data using vroom")
          wide_expression_data = vroom::vroom(tidy_expression_file) %>%
            dplyr::filter(Hugo_Symbol %in% hugo_symbols)
        } else if(engine == "fread"){
          message("Will read the data using fread")
          wide_expression_data = data.table::fread(tidy_expression_file, ) %>%
            dplyr::filter(Hugo_Symbol %in% hugo_symbols)
        } else if(engine == "grep"){
          message("Will read the data using grep")
          #lazily filter on the fly to conserve RAM (use grep without regex)
          genes_regex=paste(c("-e Hugo_Symbol",hugo_symbols),collapse = " -e ");
          grep_cmd = paste0("grep -w -F ",genes_regex," ",tidy_expression_file)
          print(grep_cmd)
          wide_expression_data = fread(cmd=grep_cmd) %>%
            dplyr::filter(Hugo_Symbol %in% hugo_symbols)
        } else {
          stop("You did not specify valid engine. Please use one of read_tsv, grep, vroom, or fread")
        }
        wide_expression_data = wide_expression_data %>%
          dplyr::select(-ensembl_gene_id) %>%
          group_by(mrna_sample_id,Hugo_Symbol) %>% #deal with non 1:1 mapping of Hugo to Ensembl
          slice_head() %>%
          as.data.frame() %>%
          pivot_wider(names_from = Hugo_Symbol, values_from = expression)
      }
      if(!missing(ensembl_gene_ids)){
        wide_expression_data = suppressMessages(read_tsv(tidy_expression_file,lazy=TRUE)) %>%
          dplyr::select(-Hugo_Symbol) %>%
          dplyr::filter(ensembl_gene_id %in% ensembl_gene_ids) %>%
          as.data.frame() %>%
          pivot_wider(names_from = ensembl_gene_id, values_from = expression)
      }
    }
  }
  
  # join the expression data to the sample ids of the metadata
  if(missing(expression_data)){
    join_with = ifelse(is.null(join_with), "NULL", join_with)
    
    if(join_with == "NULL" | join_with == "genome" | join_with == "capture"){
      these_samples_metadata = dplyr::select(these_samples_metadata, sample_id, 
                                             patient_id, biopsy_id, seq_type)
      expression_wider = left_join(these_samples_metadata, wide_expression_data,
                                   by = c("patient_id", "biopsy_id"))
      
      # add column `multi_exp` to inform whether there are more than one 
      # `mrna_sample_id` associated to a `sample_id`
      expression_wider = filter(expression_wider, !is.na(mrna_sample_id)) %>% 
        { split(.$mrna_sample_id, .$sample_id) } %>% 
        .[lengths(.) > 1] %>% 
        lapply(unique) %>% 
        .[lengths(.) > 1] %>% 
        names %>% 
        { mutate(expression_wider, multi_exp = ifelse(sample_id %in% ., 1, 0)) }
      
      ### filter out duplicated expressions based on prioritize_rows_by
      if( !missing(prioritize_rows_by) & any(expression_wider$multi_exp == 1) ){
        
        # take only duplicated gene expression rows and split them by sample_id
        multi_exp_split = dplyr::filter(expression_wider, multi_exp == 1) %>% 
          split(.$sample_id)
        
        # use the first vector of the `prioritize_rows_by` list for the filtering. 
        # if any duplication remains, use the next vector, and so on.
        multi_exp_split = lapply(multi_exp_split, function(multi_exp_split_i){
          for(based_column in names(prioritize_rows_by)){
            for( prioritize_this_value in prioritize_rows_by[[based_column]] ){
              k = multi_exp_split_i[,based_column] == prioritize_this_value
              if( any(k) ){
                multi_exp_split_i = multi_exp_split_i[k,]
                break
              } # else:
              # if there is no row with this value to prioritize, keep the rows
              # and try the lower-priority value of the next loop.
            }
          }
          # update multi_exp column
          not_duplicated = unique(multi_exp_split_i$mrna_sample_id) %>% 
            length %>% 
            {. == 1}
          if(not_duplicated){
            multi_exp_split_i$multi_exp = 0
          }
          multi_exp_split_i
        })
        
        # update expression_wider. this time duplicated rows fixed (for those sample 
        # ids that was possible)
        expression_wider = dplyr::filter(expression_wider, multi_exp == 0) %>% 
          list %>% 
          c(multi_exp_split) %>% 
          bind_rows %>% 
          arrange(sample_id, biopsy_id, patient_id, seq_type)
      }
      
    }else if(join_with == "mrna"){
      these_samples_metadata = dplyr::select(these_samples_metadata, sample_id, 
                                             patient_id, biopsy_id, seq_type)
      expression_wider = dplyr::select(wide_expression_data, -patient_id, -biopsy_id) %>% 
        left_join(these_samples_metadata, ., by = c("sample_id" = "mrna_sample_id"))
      
    }
  }else{ # expression_data is not missing
    expression_wider = wide_expression_data
  }
  
  return(expression_wider)
}
