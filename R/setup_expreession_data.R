#' @title Setup Expression Data (cBioPortal).
#'
#' @description Generate expression data based on a set of genes, format and export data for cBioPortal.
#'
#' @details This function takes a set of genes with the `these_genes` (character of vectors) parameter and returns expression data.
#' Expression data is then formatted to match the expected format for import to a cBioPortal study.
#' If no genes are provided, the function will default to all genes that are defined in the `lymphoma_genes` bundled data.
#' This function internally calls [GAMBLR.results::get_gene_expression] for returning expression data as outlined above.
#'
#' @param project_name Unique ID for your project.
#' @param clinical_file_path The path to the study specific clinical file (data_clinical_samples.txt).
#' @param these_genes Specify a set of genes (character of vectors) that you want to return expression data for. If no genes are provided, this function will resort to all lymphoma genes.
#' @param expression_df Optional argument for providing an already loaded expression matrix.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return A vector of characters with sample IDs that expression data was generated for.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' #return expression data for lymphoma genes (all samples)
#' expression_ids = setup_expression_data(out_dir = "../../")
#' }
setup_expreession_data = function(project_name = "gambl_genome",
                                  clinical_file_path = "data_clinical_samples.txt",
                                  these_genes,
                                  expression_df,
                                  out_dir){

  #get path to the clinical file holding all sample IDs
  clinical_file = data.table::fread(file = paste0(out_dir, clinical_file_path), sep = "\t", header = FALSE, skip = 5)

  #pull sample IDs
  clinical_ids = clinical_file %>%
    pull(V2)

  #default to all lymphoma genes if no specific genes are supplied for getting expression values.
  if(missing(these_genes)){
    these_genes = lymphoma_genes %>% pull(Gene)
  }

  if(missing(expression_df)){
    #get expression data
    expression_matrix = get_gene_expression(hugo_symbols = these_genes, join_with = "genome")
  }else{
    expression_matrix = expression_df
  }

  #expression metadata
  meta_expression = paste0(out_dir, "meta_expression.txt")

  meta_expression_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                   "genetic_alteration_type: MRNA_EXPRESSION\n",
                                   "datatype: CONTINUOUS\n",
                                   "stable_id: rna_seq_mrna\n",
                                   "show_profile_in_analysis_tab: false\n",
                                   "profile_name: mRNA expression (microarray)\n",
                                   "profile_description: Expression levels (Alignment microarray)\n",
                                   "data_filename: data_expressions.txt\n")

  cat(meta_expression_content, file = meta_expression)

  #subset sample IDs from the expression matrix
  expression_samples = pull(expression_matrix, sample_id) %>%
    unique()

  #intersect sample IDs from the clinical file with specified sample IDs (these_sample_ids) to
  #ensure no IDs are described in the new case list, that are not represented in the clinical file.
  ids = intersect(clinical_ids, expression_samples)

  #filter expression matrix on the sample IDs that are defined in the clinical file
  expression_matrix = expression_matrix %>% dplyr::filter(sample_id %in% ids) %>%
    dplyr::select(-capture_sample_id)

  #transform the expression matrix to match expected format
  tmp_exp <- data.frame(t(expression_matrix[]))
  names(tmp_exp) <- tmp_exp[1,]
  tmp_exp <- tmp_exp[-1,]
  expression_matrix <- tibble::rownames_to_column(tmp_exp, "Hugo_Symbol")

  #write expression matrix to file
  data_expression_full = paste0(out_dir, "data_expressions.txt")
  write_tsv(expression_matrix, data_expression_full)

  #create case list for expression data
  caselist_expression = paste0(out_dir, "case_lists/cases_expression.txt")

  tabseplist = paste(ids, collapse = "\t")

  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, "_expression"),
                   paste0("case_list_name: ", "Expression Data"),
                   paste0("case_list_description: ", "This is this case list that contains all samples that have expression data."),
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))

  cat(caselistdata, sep = "\n", file = caselist_expression)

  return(ids)
}
