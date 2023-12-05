#' @title Get Gene Expression.
#'
#' @description Get the expression for one or more genes for all GAMBL samples.
#'
#' @details Effectively get gene expression for one or multiple genes for al GAMBL samples.
#' This function can also take an already loaded expression matrix (`expression_data`)
#' to prevent the user from having to load the full expression matrix if this function needs to be run in an interactive session.
#' For examples and more info, refer to the parameter descriptions as well as vignette examples.
#' The function has argument `engine`, which accepts string "read_tsv", "grep", "vroom", and "fread". This will determine the way
#' the data is imported into R. When testing on GSC, the grep was the fastest but with a lot of variation in the run time (anywhere between 4-10 min).
#' Other engines produced similar run times (~ 7 min) on GSC, with vroom engine being the most consistent one. However, on other
#' systems (especially with fast hard drives and MacBooks for remote users) the read_tsv engine was significantly faster than the grep.
#'
#' @param metadata GAMBL metadata.
#' @param hugo_symbols One or more gene symbols. Should match the values in a maf file.
#' @param ensembl_gene_ids One or more ensembl gene IDs. Only one of hugo_symbols or ensembl_gene_ids may be used.
#' @param engine Specific way to import the data into R. Defaults to "read_tsv". Other acceptable options are "grep", "vroom", and "fread".
#' @param join_with How to restrict cases for the join. Can be one of genome, mrna or "any".
#' @param all_genes Set to TRUE to return the full expression data frame without any subsetting. Avoid this if you don't want to use tons of RAM.
#' @param expression_data Optional argument to use an already loaded expression data frame (prevent function to re-load full df from flat file or database).
#' @param from_flatfile Deprecated but left here for backwards compatibility.
#'
#' @return A data frame with gene expression.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @rawNamespace import(vroom, except = c("col_skip", "fwf_positions", "default_locale", "date_names_lang", "cols_only", "output_column", "col_character", "col_guess", "spec", "as.col_spec", "fwf_cols", "cols", "col_date", "col_datetime", "locale", "col_time", "cols_condense", "col_logical", "col_number", "col_integer", "col_factor", "fwf_widths", "date_names_langs", "problems", "date_names", "col_double", "fwf_empty"))
#' @import dplyr readr tidyr GAMBLR.data
#' @export
#'
#' @examples
#' MYC_expr = get_gene_expression(hugo_symbols = c("MYC"), join_with = "mrna")
#'
#' #Read full expression values df (no subsetting on genes)
#' full_expression_df = get_gene_expression(all_genes = TRUE,
#'                                              join_with = "genome")
#'
#' #Use loaded df (in the previous step) to get expression values for IRF4 and MYC.
#' irf4_myc_expressions = get_gene_expression(hugo_symbols = c("IRF4", "MYC"),
#'                                                all_genes = FALSE,
#'                                                join_with = "genome",
#'                                                from_flatfile = FALSE,
#'                                                expression_data = full_expression_df)
#'
get_gene_expression = function(metadata,
                               hugo_symbols,
                               ensembl_gene_ids,
                               engine = "read_tsv",
                               join_with = "mrna",
                               all_genes = FALSE,
                               expression_data,
                               from_flatfile = TRUE){
  if(!missing(hugo_symbols)){
    hugo_symbols <- as.character(hugo_symbols)
  }else if(!missing(ensembl_gene_ids)){
    ensembl_gene_ids <- as.character(ensembl_gene_ids)
  }

  database_name = GAMBLR.helpers::check_config_value(config::get("database_name"))
  if(missing(metadata)){
    if(join_with == "mrna"){
      metadata = get_gambl_metadata(seq_type_filter = "mrna", only_available = FALSE)
      metadata = metadata %>%
        dplyr::select(sample_id)

      }else if(join_with == "genome"){
      metadata = get_gambl_metadata(only_available = FALSE)
      metadata = metadata %>%
        dplyr::select(sample_id)

      }else{
      metadata = get_gambl_metadata(seq_type_filter = c("genome","mrna"), only_available = FALSE)
      metadata = metadata %>%
        dplyr::select(sample_id, biopsy_id)
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

  if(join_with == "mrna" & missing(expression_data)){
    expression_wider = dplyr::select(wide_expression_data, -biopsy_id, -genome_sample_id)
    expression_wider = left_join(metadata, expression_wider, by = c("sample_id" = "mrna_sample_id"))

    }else if(join_with == "genome" & missing(expression_data)){
    expression_wider = dplyr::select(wide_expression_data, -mrna_sample_id, -biopsy_id) %>% dplyr::filter(genome_sample_id != "NA")
    expression_wider = left_join(metadata, expression_wider, by = c("sample_id" = "genome_sample_id"))

    }else if(join_with == "any" & missing(expression_data)){
    expression_wider = dplyr::select(wide_expression_data, -mrna_sample_id, -genome_sample_id)
    expression_wider = left_join(metadata, expression_wider, by = c("biopsy_id" = "biopsy_id"))

    }else if(join_with == "mrna" & !missing(expression_data)){
      expression_wider = wide_expression_data

    }else if(join_with == "genome" & !missing(expression_data)){
      expression_wider = wide_expression_data

    }else if(join_with == "any" & !missing(expression_data)){
      expression_wider = wide_expression_data
  }
  return(expression_wider)
}
