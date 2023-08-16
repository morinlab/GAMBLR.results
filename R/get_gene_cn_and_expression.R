#' @title Get Gene CN and Expression.
#'
#' @description Get the copy number and expression for a single gene.
#'
#' @details This function works well with both Hugo Symbols and Ensembl Gene IDs.
#' It's also possible to specify more than one gene.
#'
#' @param gene_symbol One or more gene symbols. Should match the values in a maf file.
#' @param ensembl_id One or more ensembl gene IDs. Only one of hugo_symbols or ensembl_gene_ids may be used.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @return A data frame with copy number information and gene expressions.
#'
#' @import dplyr tibble
#' @export
#'
#' @examples
#' MYC_cn_expression = get_gene_cn_and_expression("MYC")
#'
get_gene_cn_and_expression = function(gene_symbol,
                                      ensembl_id,
                                      this_seq_type = "genome"){

    if(!missing(gene_symbol)){
      this_row = GAMBLR.data::grch37_gene_coordinates %>%
        dplyr::filter(hugo_symbol == gene_symbol)

      this_region = paste0(this_row$chromosome, ":", this_row$start, "-", this_row$end)
      gene_name = gene_symbol
      }

    else{
      this_row = GAMBLR.data::grch37_gene_coordinates %>%
        dplyr::filter(ensembl_gene_id == ensembl_id)

      this_region = paste0(this_row$chromosome, ":", this_row$start, "-",this_row$end)
      gene_name = ensembl_id
      gene_symbol = pull(this_row, hugo_symbol)
    }
  gene_cn = get_cn_states(
      regions_list = c(this_region),
      region_names = c(gene_name),
      this_seq_type = this_seq_type) %>%
    as.data.frame()

  colnames(gene_cn)[1] = paste(colnames(gene_cn)[1], "CN", sep = "_")
  gene_cn = gene_cn %>%
    rownames_to_column("sample_id")

  gene_exp = get_gene_expression(hugo_symbols = c(gene_symbol), join_with = "genome")
  exp_copy = left_join(gene_cn, gene_exp, by = "sample_id")
  all_meta = get_gambl_metadata()
  exp_copy_meta = left_join(all_meta, exp_copy, by = "sample_id")
  return(exp_copy_meta)
}
