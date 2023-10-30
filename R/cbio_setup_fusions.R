#' @title Setup Fusions (cBioPortal).
#'
#' @description Annotate SVs and create the input for fusions to be displayed in cBioPortal instance.
#'
#' @details This function calls [GAMBLR.results::get_combined_sv] and runs [GAMBLR.utils::annotate_sv] on the returned data frame.
#' Should be run as the next step after running [GAMBLR.results::setup_study]. Note that the parameters called with this function
#' has to match the parameter arguments of [GAMBLR.results::setup_study], i.e if `short_name` is for [GAMBLR.results::setup_study] is "GAMBL",
#' then the `short_name` in [GAMBLR.results::setup_fusions] also has to be "GAMBL", etc.
#'
#' @param short_name A concise name for your portal project.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param gambl_maf maf origin.
#' @param gambl_icgc_maf ICGC maf origin.
#' @param description A verbose description of your data set.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return A vector of sample_id for the patients that have been included.
#'
#' @import dplyr readr tidyr GAMBLR.utils
#' @export
#'
#' @examples
#' \dontrun{
#' fusion_ids = cbio_setup_fusions(out_dir = "GAMBLR/cBioPortal/instance01/")
#' }
cbio_setup_fusions = function(short_name = "GAMBL",
                              human_friendly_name = "GAMBL data",
                              project_name = "gambl_genome",
                              description = "GAMBL data from genome",
                              gambl_maf = "maf_slms3_hg19",
                              gambl_icgc_maf = "maf_slms3_hg19_icgc",
                              out_dir){


  #create necessary files
  #meta fusions
  meta_fusions = paste0(out_dir, "meta_fusions.txt")

  meta_fusion_content = paste0("cancer_study_identifier: ", project_name, "\n",
                               "genetic_alteration_type: FUSION\n",
                               "datatype: FUSION\n",
                               "stable_id: fusion\n",
                               "show_profile_in_analysis_tab: true\n",
                               "profile_name: Fusions\n",
                               "profile_description: Fusion data\n",
                               "data_filename: data_fusions.txt\n")

  cat(meta_fusion_content, file = meta_fusions)

  #get SV breakpoints and annotate them
  unannotated_sv = get_combined_sv()

  annotated_sv = GAMBLR.utils::annotate_sv(unannotated_sv) %>%
    dplyr::filter(!is.na(partner)) %>%
    as.data.frame()

  fusion_samples = pull(annotated_sv, tumour_sample_id) %>%
    unique()

  #deal with any cases not in metadata
  fusions_df =  data.frame(Hugo_Symbol = annotated_sv$gene,
                           Center = "BCGSC",
                           Tumor_Sample_Barcode = annotated_sv$tumour_sample_id,
                           Fusion = c(pull(unite(annotated_sv, fusion, partner, gene, sep = "-"), fusion)),
                           DNA_support = "yes",
                           RNA_support = "no",
                           Method = "SVAR",
                           Frame = "in-frame")

  fusions_df = distinct(fusions_df, Tumor_Sample_Barcode, Fusion, .keep_all = TRUE)

  #determine what table to query and what restrictions to use for the MAF data

  # TO DO: Fix this code to work with the indexed MAF file using get_ssm_by_region instead of by gene
  #nfkbiz_entrez = 64332
  #nfkbiz_utr_ssm = get_ssm_by_gene(gene_symbol = "NFKBIZ") %>%
  #  dplyr::filter(Variant_Classification == "3'UTR") %>%
  #  pull(Tumor_Sample_Barcode) %>%
  #  unique()

  #nfkbiz.mut.df = data.frame(Hugo_Symbol = "NFKBIZ",
  #                           Entrez_Gene_Id = nfkbiz_entrez,
  #                           Center = "BCGSC",
  #                           Tumor_Sample_Barcode = nfkbiz_utr_ssm,
  #                           Fusion = "NFKBIZ-UTR",
  #                           DNA_support = "yes",
  #                           RNA_support = "no",
  #                           Method = "SLMS-3",
  #                           Frame = "in-frame")

  #get any SV breakpoints that are in the 3'UTR of NFKBIZ
  #nfkbiz_utr_region = "chr3:101,578,185-101,579,902"
  data_fusions = paste0(out_dir, "data_fusions.txt")
  #TODO: FIX! this is also broken
  #nfkbiz.svs = get_combined_sv(region = nfkbiz_utr_region) %>%
  #  pull(tumour_sample_id) %>%
  #  unique()

  #nfkbiz.sv.df = data.frame(Hugo_Symbol = "NFKBIZ",
  #                          Entrez_Gene_Id = nfkbiz_entrez,
  #                          Center = "BCGSC",
  #                          Tumor_Sample_Barcode = nfkbiz.svs,
  #                          Fusion = "NFKBIZ-SV",
  #                          DNA_support = "yes",
  #                          RNA_support = "no",
  #                          Method = "Manta",
  #                          Frame = "in-frame")

  #all_fusions = rbind(fusions_df, nfkbiz.sv.df, nfkbiz.mut.df)
  all_fusions = fusions_df
  fusion.cases = as.character(unique(all_fusions$Tumor_Sample_Barcode))
  write_tsv(all_fusions, data_fusions)

  #create necessary files
  #create caselist fusions
  caselist_fusion = paste0(out_dir, "case_lists/cases_fusion.txt")

  tabseplist = paste(fusion.cases, collapse = "\t")

  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, "_fusions"), "case_list_name: Samples with fusions.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))

  cat(caselistdata, sep = "\n", file = caselist_fusion)
  return(fusion.cases)
}
