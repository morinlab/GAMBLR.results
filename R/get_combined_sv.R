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
#' [GAMBLR.results::get_manta_sv], [GAMBLR.results::get_manta_sv_by_sample], [GAMBLR.results::get_manta_sv_by_samples]
#'
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0. (default: 0)
#' @param these_sample_ids A character vector of tumour sample IDs you wish to retrieve SVs for.
#' @param these_samples_metadata a GAMBL metadata frame containing only the samples you want returned
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses). Default is FALSE.
#' @param projection The projection genome build. Default is "grch37".
#' @param oncogenes A character vector of genes commonly involved in translocations.
#' You may want to include one or more e.g.: c("CCND1", "BCL2","MYC")
#' NOTE: If you are looking for SV affecting oncogenes you are more likely
#' going to want to pass the full output to [GAMBLR.utils::annotate_sv]
#' @param region Optional, region formatted like chrX:1234-5678 (chromosome can be prefixed or not prefixed) instead of specifying chromosome, start and end separately.
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs.
#'
#' @import config dplyr readr stringr GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#'   all_sv_bed = get_combined_sv()
#'   annotated_onco_sv = GAMBLR.utils::annotate_sv(all_sv_bed) %>% 
#'         dplyr::filter(gene %in% c("BCL2", "BCL6", "MYC", "CCND1"),!is.na(partner))
#' }
get_combined_sv = function(min_vaf = 0,
                           these_sample_ids,
                           these_samples_metadata,
                           with_chr_prefix = FALSE,
                           projection = "grch37",
                           oncogenes,
                           region){
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }
  
  base_path = GAMBLR.helpers::check_config_value(config::get("project_base"))
  sv_file = GAMBLR.helpers::check_config_value(config::get()$results_flatfiles$sv_combined$icgc_dart)
  if(projection == "hg38"){
    sv_file = str_replace(sv_file, "--grch37", "--hg38")
  }
  sv_file = paste0(base_path, sv_file)
  permissions = file.access(sv_file, 4)
  if(permissions == - 1){
    sv_file = GAMBLR.helpers::check_config_value(config::get()$results_flatfiles$sv_combined$gambl)
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
  if(!missing(these_samples_metadata)){
    all_sv = all_sv %>%
      dplyr::filter(tumour_sample_id %in% these_samples_metadata$sample_id)
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
  #deal with chr prefixes based on the selected projection (if return is to be subset to regions...)
  if(!missing(region)){
    if(projection == "grch37"){
      if(grepl("chr", chromosome)){
        chromosome = gsub("chr", "", chromosome)
      }
    }else if(projection == "hg38"){
      if(!grepl("chr", chromosome)){
        chromosome = paste0("chr", chromosome)
      }
    }
    all_sv = all_sv %>%
      dplyr::filter((CHROM_A == chromosome & START_A >= qstart & START_A <= qend) | (CHROM_B == chromosome & START_B >= qstart & START_B <= qend))
  }
  
  all_sv = all_sv %>%
      mutate(FILTER = "PASS") #i.e all variants returned with get_combined_sv() all have PASS in the FILTER column.

  return(all_sv)
}
