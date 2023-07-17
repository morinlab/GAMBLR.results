#' @title Get Combined SV.
#'
#' @description Retrieve Combined Manta and GRIDSS-derived SVs from a flatfile and filter.
#'
#' @details The bedpe files used as input to this function were pre-filtered for a minimum VAF of 0.05, and SVs affecting.
#' common translocation regions (BCL2, BCL6, MYC, CCND1) were whitelisted (e.g. no VAF filter applied).
#' Therefore if you wish to post-filter the SVs we recommend doing so carefully after loading this data frame.
#' Further, the input bedpe file is annotated with oncogenes and superenhancers from naive and germinal centre B-cells.
#' You can subset to events affecting certain loci using the "oncogenes" argument.
#' Is this function not what you are looking for? Try one of the following, similar, functions;
#' [GAMBLR::get_manta_sv], [GAMBLR::get_manta_sv_by_sample], [GAMBLR::get_manta_sv_by_samples]
#'
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0. (default: 0)
#' @param these_sample_ids A character vector of tumour sample IDs you wish to retrieve SVs for.
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses). Default is FALSE.
#' @param projection The projection genome build. Default is "grch37".
#' @param oncogenes A character vector of genes commonly involved in translocations. Possible values: CCND1, CIITA, SOCS1, BCL2, RFTN1, BCL6, MYC, PAX5.
#'
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs.
#'
#' @import config dplyr readr stringr
#' @export
#'
#' @examples
#' get_combined_sv(oncogenes = c("MYC", "BCL2", "BCL6"))
#'
get_combined_sv = function(min_vaf = 0,
                           these_sample_ids,
                           with_chr_prefix = FALSE,
                           projection = "grch37",
                           oncogenes){

  base_path = check_config_value(config::get("project_base"))
  sv_file = check_config_value(config::get()$results_flatfiles$sv_combined$icgc_dart)
  if(projection == "hg38"){
    sv_file = str_replace(sv_file, "--grch37", "--hg38")
  }
  sv_file = paste0(base_path, sv_file)
  permissions = file.access(sv_file, 4)
  if(permissions == - 1){
    sv_file = check_config_value(config::get()$results_flatfiles$sv_combined$gambl)
    sv_file = paste0(base_path, sv_file)
  }

  #check for missingness
  if(!file.exists(sv_file)){
    print(paste("missing: ", sv_file))
    message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
    message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
  }

  all_sv = suppressMessages(read_tsv(sv_file, col_types = "cnncnncnccccnnccncn")) %>%
    dplyr::rename(c("VAF_tumour" = "VAF")) %>%
    dplyr::filter(VAF_tumour >= min_vaf)

  if(!missing(these_sample_ids)){
    all_sv = all_sv %>%
      dplyr::filter(tumour_sample_id %in% these_sample_ids)
  }

  if(!missing(oncogenes)){
    all_sv = all_sv %>%
      dplyr::filter(ANNOTATION_A %in% oncogenes | ANNOTATION_B %in% oncogenes)
  }

  if(with_chr_prefix){
    #add chr prefix only if it's missing
    all_sv = all_sv %>%
      dplyr::mutate(CHROM_A = case_when(str_detect(CHROM_A, "chr") ~ CHROM_A, TRUE ~ paste0("chr", CHROM_A)))

    all_sv = all_sv %>%
      dplyr::mutate(CHROM_B = case_when(str_detect(CHROM_B, "chr") ~ CHROM_B, TRUE ~ paste0("chr", CHROM_B)))
  }

  all_sv = all_sv %>%
      mutate(FILTER = "PASS") #i.e all variants returned with get_combined_sv() all have PASS in the FILTER column.

  return(all_sv)
}
