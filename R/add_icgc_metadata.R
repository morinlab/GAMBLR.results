#' @title Add ICGC metadata.
#'
#' @description Layer on ICGC metadata from a supplemental table to fill in missing COO.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::get_gambl_metadata], not meant for out-of-package usage.
#'
#' @param incoming_metadata A metadata table (probably output from `get_gambl_metadata`).
#'
#' @return Metadata with layered information (ICGC).
#'
#' @import dplyr readr stringr
#'
#' @noRd
#'
#' @examples
#' my_meta = get_gambl_metadata()
#' icgc_metadata = add_icgc_metadata(incoming_metadata = my_meta)
#'
#' @export
add_icgc_metadata = function(incoming_metadata){
  repo_base = check_config_value(config::get("repo_base"))
  icgc_publ_file = paste0(repo_base,"data/metadata/raw_metadata/MALY_DE_tableS1.csv")
  icgc_publ = suppressMessages(suppressWarnings(read_csv(icgc_publ_file)))
  icgc_publ = icgc_publ[,c(1:20)]
  #fix commas as decimals
  icgc_publ = mutate(icgc_publ, purity = str_replace(purity, ",", "."))
  icgc_publ = mutate(icgc_publ, sex = str_to_upper(sex))

  icgc_raw_path = paste0(repo_base,"data/metadata/raw_metadata/ICGC_MALY_seq_md.tsv")

  #check for missingness
  if(!file.exists(icgc_raw_path)){
    print(paste("missing: ", icgc_raw_path))
    message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
  }

  icgc_raw = suppressMessages(read_tsv(icgc_raw_path))

  icgc_raw = icgc_raw %>%
    dplyr::select(-compression, -bam_available, -read_length, -time_point, -unix_group, -ffpe_or_frozen, -link_name)  %>%
    dplyr::filter(tissue_status %in% c("tumor", "tumour"))

  icgc_all = left_join(icgc_raw, icgc_publ,by = "ICGC_ID") %>%
    dplyr::select(-tissue_status, -seq_type, -protocol, -seq_source_type, -data_path, -genome_build, -RNA_available) %>%
    dplyr::select(sample_id, ICGC_ID, pathology.x, pathology.y, COO, molecular_BL, MYC_sv, BCL2_sv, BCL6_sv) %>%
    dplyr::rename("ICGC_MYC_sv" = "MYC_sv") %>%
    dplyr::rename("ICGC_BCL2_sv" = "BCL2_sv") %>%
    dplyr::rename("ICGC_BCL6_sv" = "BCL6_sv") %>%
    dplyr::rename("detailed_pathology" = "pathology.x") %>%
    dplyr::rename("ICGC_PATH" = "pathology.y")

  #join with all metadata to fill in blanks
  #all_meta=get_gambl_metadata()
  rejoined = left_join(incoming_metadata, icgc_all,by = "sample_id") %>%
    mutate(COO_final = case_when(!is.na(COO_consensus) ~ COO_consensus, COO != "n.a." & COO != "TypeIII" ~ COO, TRUE ~ "NA")) %>%
    dplyr::select(-COO)
  return(rejoined)
}
