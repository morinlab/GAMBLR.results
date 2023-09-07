#' @title Tidy gene Expression.
#'
#' @description Read a full expression matrix.
#'
#' @details Read a full expression matrix and subset to samples in GAMBL that have metadata (remove duplicates with consistent preferences).
#' The user can also specify if they want the data frame returned into their R session, or if the data frame should be written to file (default).
#'
#' @param return_df Boolean parameter to return the dataframe, default is FALSE (i.e writing results to file).
#'
#' @import dplyr readr stringr tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' #return data frame with gene expression to R
#' gene_expression = tidy_gene_expression(return_df = TRUE)
#' }
#' 
tidy_gene_expression = function(return_df = FALSE){

  #read in the full matrix
  ex_matrix_file = GAMBLR.helpers::check_config_value(config::get("results_merged")$ex_matrix_file)
  tidy_expression_file = GAMBLR.helpers::check_config_value(config::get("results_merged")$tidy_expression_file)
  print("Loading and tidying the full matrix file...")

  ex_tidy = suppressMessages(read_tsv(ex_matrix_file)) %>%
    dplyr::select(-gene_id) %>%
    dplyr::rename("Hugo_Symbol" = "hgnc_symbol") %>%
    pivot_longer(-c(Hugo_Symbol, ensembl_gene_id), names_to = "sample_id", values_to = "expression")

  all_samples = distinct(ex_tidy, sample_id) %>%
    pull(sample_id)
  #retrieve the full list of sample_id for RNA-seq libraries that have data in this matrix

  #pull the full metadata for all RNA-seq samples in GAMBL
  #under the hood this is a join of the sample and biopsy tables but subset for RNA-seq
  print("Loading the mrna metadata...")
  rna_meta = get_gambl_metadata(seq_type_filter = "mrna", tissue_status_filter = c("tumour", "normal"), only_available = FALSE)
  #subset to just the ones in the matrix and keep only the relevant rows
  rna_meta_existing = rna_meta %>%
    dplyr::filter(sample_id %in% all_samples) %>%
    dplyr::select(sample_id,patient_id, biopsy_id, protocol, cohort, tissue_status) %>%
    mutate(biopsy_id = ifelse(tissue_status == "normal", sample_id, biopsy_id))

  print("Selecting one library per biopsy...")
  selected_libraries = rna_meta_existing %>%
    group_by(biopsy_id) %>%
    # Take the biopsy_id with the longest string length (e.g PolyA vs. Ribodepletion)
    slice_max(str_length(protocol), n = 1, with_ties = FALSE) %>%
    ungroup()
  #this brings the total from 1423 to 1322

  ex_tidy = ex_tidy %>%
    dplyr::rename(mrna_sample_id = sample_id)

  ex_tidy = ex_tidy %>%
    dplyr::filter(mrna_sample_id %in% selected_libraries$sample_id)

  rna_meta = rna_meta %>%
    dplyr::select(sample_id, biopsy_id)

  ex_tidy_final = left_join(ex_tidy, rna_meta, by = c("mrna_sample_id" = "sample_id"))
  #this still has the mrna sample ID. Need to add the genome sample_id from this biopsy (where available)
  print("Loading the genome metadata...")
  genome_meta = get_gambl_metadata() %>%
    dplyr::select(biopsy_id, sample_id, patient_id, ffpe_or_frozen) %>%
    dplyr::rename("genome_sample_id" = "sample_id")
  #this gets the metadata in the same format but restricted to genome samples

  #join to genome metadata based on biopsy_id (should be the same for RNA-seq and tumour genomes)
  ex_tidy_genome = left_join(ex_tidy_final, genome_meta, by = "biopsy_id")

  ex_tidy_genome = dplyr::select(ex_tidy_genome, ensembl_gene_id, Hugo_Symbol, mrna_sample_id, expression, biopsy_id, genome_sample_id)

  #write the data back out for use by others and loading into the database.
  if(return_df){
    return(ex_tidy_genome)
  }else{
    print("Writing the tidy data to file...")
    write_tsv(ex_tidy_genome, file = tidy_expression_file)
  }
  print("Done!")
}
