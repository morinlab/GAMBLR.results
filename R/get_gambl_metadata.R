#' @title Get GAMBL metadata.
#'
#' @description Return metadata for a selection of samples.
#'
#' @details This function returns metadata for GAMBL samples. Options for subset and filter the returned data are available.
#' For more information on how to use this function with different filtering criteria, refer to the parameter descriptions,
#' examples and vignettes. Embargoed cases (current options: 'BLGSP-study', 'FL-study', 'DLBCL-study', 'FL-DLBCL-study', 'FL-DLBCL-all', 'DLBCL-unembargoed', 'BL-DLBCL-manuscript', 'MCL','MCL-CLL')
#'
#' @param seq_type_filter Filtering criteria (default: all genomes).
#' @param tissue_status_filter Filtering criteria (default: only tumour genomes, can be "mrna" or "any" for the superset of cases).
#' @param case_set Optional short name for a pre-defined set of cases avoiding any embargoed cases (current options: 'BLGSP-study', 'FL-study', 'DLBCL-study', 'FL-DLBCL-study', 'FL-DLBCL-all', 'DLBCL-unembargoed', 'BL-DLBCL-manuscript', 'MCL','MCL-CLL').
#' @param remove_benchmarking By default the FFPE benchmarking duplicate samples will be dropped.
#' @param sample_flatfile Optionally provide the full path to a sample table to use instead of the default.
#' @param biopsy_flatfile Optionally provide the full path to a biopsy table to use instead of the default.
#' @param with_outcomes Optionally join to gambl outcome data.
#' @param only_available If TRUE, will remove samples with FALSE or NA in the bam_available column (default: TRUE).
#' @param from_flatfile New default is to use the metadata in the flat-files from your clone of the repo. Can be overridden to use the database.
#' @param seq_type_priority For duplicate sample_id with different seq_type available, the metadata will prioritize this seq_type and drop the others.
#'
#' @return A data frame with metadata for each biopsy in GAMBL
#'
#' @import config dplyr tidyr readr RMariaDB DBI
#' @export
#'
#' @examples
#' #basic usage
#' my_metadata = get_gambl_metadata()
#'
#' #use pre-defined custom sample sets
#' only_blgsp_metadata = get_gambl_metadata(case_set = "BLGSP-study")
#'
#' #override default filters and request metadata for samples other than tumour genomes,
#' #e.g. also get the normals
#' only_normal_metadata = get_gambl_metadata(tissue_status_filter = c('tumour','normal'))
#'
#' non_duplicated_genome_and_capture = get_gambl_metadata(seq_type_filter = c('genome', 'capture'),
#'                                                        seq_type_priority = "genome")
#'
get_gambl_metadata = function(seq_type_filter = "genome",
                              tissue_status_filter = "tumour",
                              case_set,
                              remove_benchmarking = TRUE,
                              with_outcomes = TRUE,
                              from_flatfile = TRUE,
                              sample_flatfile,
                              biopsy_flatfile,
                              only_available = TRUE,
                              seq_type_priority = "genome"){

  check_remote_configuration()
  #this needs to be in any function that reads files from the bundled GAMBL outputs synced by Snakemake
  outcome_table = get_gambl_outcomes(from_flatfile = from_flatfile)

  if(from_flatfile){
    base = config::get("repo_base")
    if(missing(sample_flatfile)){
      sample_flatfile = paste0(base, config::get("table_flatfiles")$samples)
    }
    if(missing(biopsy_flatfile)){
      biopsy_flatfile = paste0(base, config::get("table_flatfiles")$biopsies)
    }

    #check for missingness
    if(!file.exists(sample_flatfile)){
      print(paste("missing: ", sample_flatfile))
      message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
    }

    if(!file.exists(biopsy_flatfile)){
      print(paste("missing: ", biopsy_flatfile))
      message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
    }

    sample_meta = suppressMessages(read_tsv(sample_flatfile, guess_max = 100000))
    biopsy_meta = suppressMessages(read_tsv(biopsy_flatfile, guess_max = 100000))

  }else{
    db = check_config_value(config::get("database_name"))
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
    sample_meta = dplyr::tbl(con, "sample_metadata") %>%
      as.data.frame()

    biopsy_meta = dplyr::tbl(con, "biopsy_metadata") %>%
      as.data.frame()

    DBI::dbDisconnect(con)
  }

  # Conditionally remove samples without bam_available == TRUE
  if(only_available == TRUE){
    sample_meta = dplyr::filter(sample_meta, bam_available %in% c(1, "TRUE"))
  }
  sample_meta_normal_genomes =  sample_meta %>%
    dplyr::filter(seq_type %in% seq_type_filter & tissue_status == "normal") %>%
    dplyr::select(patient_id, sample_id, seq_type, genome_build) %>% as.data.frame() %>%
    #dplyr::select(patient_id, sample_id, seq_type) %>% as.data.frame() %>%
    dplyr::rename("normal_sample_id" = "sample_id")
#print(head(sample_meta_normal_genomes))
  sample_meta = sample_meta %>%
    dplyr::filter(seq_type %in% seq_type_filter & tissue_status %in% tissue_status_filter) %>%
    dplyr::select(-sex)

  #if only normals were requested, just return what we have because there is nothing else to join
  if(tissue_status_filter == "normal"){
    return(sample_meta)
  }

  #if we only care about genomes, we can drop/filter anything that isn't a tumour genome
  #The key for joining this table to the mutation information is to use sample_id. Think of this as equivalent to a library_id. It will differ depending on what assay was done to the sample.
  biopsy_meta = biopsy_meta %>%
    dplyr::select(-patient_id) %>%
    dplyr::select(-pathology) %>%
    dplyr::select(-time_point) %>%
    dplyr::select(-EBV_status_inf) #drop duplicated columns

  all_meta = dplyr::left_join(sample_meta, biopsy_meta, by = "biopsy_id") %>%
    as.data.frame()

  all_meta = all_meta %>%
    mutate(bcl2_ba = ifelse(bcl2_ba == "POS_BCC", "POS", bcl2_ba))

  if(!"mrna" %in% seq_type_filter & length(tissue_status_filter) == 1 & tissue_status_filter[1] == "tumour"){
    #join back the matched normal genome
    #all_meta = left_join(all_meta, sample_meta_normal_genomes, by=c("patient_id", "seq_type"))
    all_meta = left_join(all_meta, sample_meta_normal_genomes, by=c("patient_id", "seq_type","genome_build"))
    all_meta = all_meta %>%
      mutate(pairing_status = case_when(is.na(normal_sample_id)~"unmatched", TRUE~"matched"))
  }
  #all_meta[all_meta$pathology == "B-cell unclassified","pathology"] = "HGBL"  #TODO fix this in the metadata
  if(remove_benchmarking){
    all_meta = all_meta %>%
      dplyr::filter(cohort != "FFPE_Benchmarking")
  }
  if("any" %in% seq_type_filter){
    #this option may not work. To be deprecated, probably

   all_meta = all_meta %>%
    arrange(seq_type) %>%
    group_by(biopsy_id) %>%
    slice_head()
  }
  all_meta = add_icgc_metadata(all_meta) %>%
    mutate(consensus_pathology = case_when(ICGC_PATH == "FL-DLBCL" ~ "COM",
                                           ICGC_PATH == "DH-BL" ~ pathology,
                                           ICGC_PATH == "FL" | ICGC_PATH == "DLBCL" ~ ICGC_PATH,
                                           pathology == "COMFL" ~ "COM",
                                           TRUE ~ pathology))

  all_meta = unique(all_meta) #something in the ICGC code is causing this. Need to figure out what #should this be posted as an issue on Github?
  if(!missing(case_set)){
    # This functionality is meant to eventually replace the hard-coded case sets
    case_set_path = check_config_value(config::get("sample_sets")$default)
    full_case_set_path =  paste0(check_config_value(config::get("repo_base")), case_set_path)
    if (file.exists(full_case_set_path)) {
      full_case_set = suppressMessages(read_tsv(full_case_set_path))
    } else {
      message(paste("Warning: case set is requested, but the case set file", full_case_set_path, "is not found."))
      message("Defaulting to pre-defined case sets")
    }

    # pre-defined case sets
    if(case_set == "MCL"){
      all_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("MCL"))

    }else if(case_set == "MCL-CLL"){
      all_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("MCL", "CLL")) %>%
        dplyr::filter(cohort != "CLL_LSARP_Trios")
    }else if(case_set == "tFL-study"){
      #update all DLBCLs in this file to indicate they're transformations
      transformed_manual <- paste0(
          base,
          "data/metadata/raw_metadata/gambl_tFL_manual.tsv"
      )
      transformed_manual <- suppressMessages(
          read_tsv(
              transformed_manual
          )
      )

      all_meta = left_join(all_meta, transformed_manual)
      fl_meta_kridel = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        group_by(patient_id) %>%
        mutate(FL = sum(pathology == "FL"), DLBCL = sum(pathology %in% c("COM", "DLBCL", "COMFL"))) %>%
        mutate(transformed = ifelse(FL > 0 & DLBCL > 0, TRUE, FALSE))  %>%
        mutate(analysis_cohort = case_when(consensus_pathology == "FL" & transformed == TRUE ~ "pre-HT",
                                           consensus_pathology == "DLBCL" & transformed == TRUE ~ "ignore",
                                           TRUE ~ "no-HT")) %>%
        dplyr::filter(cohort == "FL_Kridel") %>%
        dplyr::filter((analysis_cohort == "no-HT" & time_point == "A")|(analysis_cohort == "pre-HT")) %>%
        dplyr::select(-transformed, -FL, -DLBCL)

      dlbcl_meta_kridel = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        group_by(patient_id) %>%
        mutate(FL = sum(pathology == "FL"), DLBCL = sum(pathology %in% c("COM", "DLBCL", "COMFL"))) %>%
        mutate(transformed = ifelse(FL > 0 & DLBCL > 0, TRUE, FALSE))  %>%
        mutate(analysis_cohort = case_when(consensus_pathology == "FL" & transformed == TRUE ~ "post-HT",
                                           consensus_pathology == "DLBCL" & transformed == TRUE ~ "ignore",
                                           TRUE ~ "post-HT")) %>%
        dplyr::filter((analysis_cohort == "post-HT" & time_point == "B"))

      fl_meta_other = all_meta %>% dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        dplyr::filter(cohort != "FL_Kridel") %>%
        dplyr::filter((consensus_pathology %in% c("FL", "COM"))) %>% mutate(analysis_cohort = consensus_pathology)

      gambl_transformations <- paste0(
          base,
          "data/metadata/raw_metadata/gambl_transformation.txt"
      )
      gambl_transformations <- suppressMessages(
              read_delim(
                  gambl_transformations,
                  delim = " "
              )
          ) %>%
          dplyr::filter(code_transf == 1) %>%
          group_by(res_id) %>%
          slice_head()

      fl_transformation_meta = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/shared/gambl_fl_transformed.tsv"))
      transformed_cases = pull(gambl_transformations, res_id)
      fl_meta_other[which(fl_meta_other$patient_id %in% transformed_cases), "analysis_cohort"] = "pre-HT"
      fl_meta_other = mutate(fl_meta_other, analysis_cohort = ifelse(analysis_cohort == "FL", "no-HT", analysis_cohort))

      #Finally over-ride analysis cohort with the outcome of clinical review
      dlbcl_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP", "FL_Kridel", "FFPE_Benchmarking")) %>%
        dplyr::filter(consensus_pathology == "DLBCL") %>%
        mutate(analysis_cohort = "denovo-DLBCL")

      all_meta  = bind_rows(dlbcl_meta, dlbcl_meta_kridel, fl_meta_kridel, fl_meta_other) %>%
        unique()

      curated <- paste0(
          base,
          "data/metadata/raw_metadata/clin_review_fl.tsv"
      )
      curated <- suppressMessages(
          read_tsv(
              curated
          )
      )

      all_meta = left_join(all_meta, curated) %>%
        mutate(analysis_cohort = ifelse(is.na(clinical_review), analysis_cohort, clinical_review))

      all_meta[which(all_meta$is_tFL == 1), "analysis_cohort"] = "post-HT"
    } else if(case_set == "FL-DLBCL-study"){

      #get FL cases and DLBCL cases not in special/embargoed cohorts
      fl_meta_kridel = all_meta %>% dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        group_by(patient_id) %>%
        mutate(FL = sum(pathology == "FL"), DLBCL = sum(pathology %in% c("COM", "DLBCL", "COMFL"))) %>%
        mutate(transformed = ifelse(FL > 0 & DLBCL > 0, TRUE, FALSE)) %>%
        mutate(analysis_cohort=case_when(consensus_pathology == "FL" & transformed == TRUE ~ "tFL", consensus_pathology == "DLBCL" & transformed == TRUE ~ "ignore", TRUE ~ "FL")) %>%
        dplyr::filter(cohort == "FL_Kridel") %>%
        dplyr::filter((analysis_cohort == "FL" & time_point == "A")|(analysis_cohort == "tFL")) %>%
        dplyr::select(-transformed, -FL, -DLBCL)

      fl_meta_other = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        dplyr::filter(cohort != "FL_Kridel") %>%
        dplyr::filter((consensus_pathology %in% c("FL", "COM"))) %>% mutate(analysis_cohort = consensus_pathology)

      fl_transformation_meta = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/shared/gambl_fl_transformed.tsv"))
      transformed_cases = fl_transformation_meta %>%
        dplyr::filter(!is.na(PATHa.tr)) %>%
        pull(patient_id)

      fl_meta_other[which(fl_meta_other$patient_id %in% transformed_cases), "analysis_cohort"] = "tFL"

      dlbcl_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP", "FL_Kridel", "FFPE_Benchmarking")) %>%
        dplyr::filter(consensus_pathology == "DLBCL" & COO_final == "GCB") %>%
        mutate(analysis_cohort = "DLBCL")

      all_meta = bind_rows(dlbcl_meta, fl_meta_kridel, fl_meta_other) %>%
      unique()

    }else if(case_set == "FL-study"){

      #get FL cases and DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios")) %>%
        group_by(patient_id) %>%
        arrange(patient_id, pathology)  %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::filter(pathology == "FL")

    }else if(case_set == "DLBCL-study"){

      #get FL cases and DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios")) %>%
        group_by(patient_id) %>%
        arrange(patient_id, pathology)  %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::filter(consensus_pathology %in% c("DLBCL", "FL"))

    }else if(case_set == "DLBCL-unembargoed"){
      #get DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>%
      dplyr::filter(consensus_pathology %in% c("DLBCL", "COM")) %>%
      dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP"))

    }else if(case_set == "BLGSP-study"){

      #get BL cases minus duplicates (i.e. drop benchmarking cases)
      all_meta = all_meta %>%
        dplyr::filter(cohort %in% c("BL_Adult", "BL_cell_lines", "BL_ICGC", "BLGSP_Bcell_UNC", "BL_Pediatric") | (sample_id == "06-29223T"))

    }else if(case_set == "BL-DLBCL-manuscript"){
      set_file = paste0(base, "data/metadata/BLGSP--DLBCL-case-set.tsv")
      adult_bl_manuscript_samples = read_tsv(set_file) %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
      dplyr::filter(sample_id %in% adult_bl_manuscript_samples)

    }else if(case_set == "BL-DLBCL-manuscript-HTMCP"){
      set_file = paste0(base, "data/metadata/BLGSP--DLBCL-case-set.tsv")
      adult_bl_manuscript_samples = read_tsv(set_file) %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
      dplyr::filter(sample_id %in% adult_bl_manuscript_samples | cohort == "DLBCL_HTMCP")

    }else if(case_set == "FL-DLBCL-all"){
      set_file = paste0(base, "data/metadata/FL--DLBCL--all-case-set.tsv")
      fl_dlbcl_all_samples = read_tsv(set_file) %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
      dplyr::filter(sample_id %in% fl_dlbcl_all_samples)

    }else if(case_set == "GAMBL-all"){

      #get all GAMBL but remove FFPE benchmarking cases and ctDNA
      all_meta = all_meta %>%
      dplyr::filter(!cohort %in% c("FFPE_Benchmarking", "DLBCL_ctDNA"))
    }else if(case_set %in% colnames(full_case_set)) {
      # ensure consistent column naming
      full_case_set =
        full_case_set %>% rename_at(vars(matches(
          "sample_id", ignore.case = TRUE
        )),
        ~ "Tumor_Sample_Barcode")

      # get case set as defined in the file
      this_subset_samples =
        full_case_set %>%
        dplyr::filter(!!sym(case_set) == 1) %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
        dplyr::filter(sample_id %in% this_subset_samples)

    }else{
      message(paste("case set", case_set, "not available"))
      return()
    }
  }

  all_meta = GAMBLR::tidy_lymphgen(all_meta,
              lymphgen_column_in = "lymphgen_cnv_noA53",
              lymphgen_column_out = "lymphgen",
              relevel=TRUE)

  #all_meta = GAMBLR::collate_lymphgen(all_meta, verbose=FALSE)

  # "catchall" pathology for those that need review
  all_meta = all_meta %>%
    mutate(pathology = ifelse(nchar(pathology) > 15, "OTHER", pathology))

  all_meta = mutate(all_meta, Tumor_Sample_Barcode = sample_id) #duplicate for convenience
  all_meta = all_meta %>%
    dplyr::mutate(consensus_coo_dhitsig = case_when(pathology != "DLBCL" ~ pathology,
                                                    COO_consensus == "ABC" ~ COO_consensus,
                                                    DLBCL90_dhitsig_call == "POS" ~ "DHITsigPos",
                                                    DLBCL90_dhitsig_call == "NEG" ~ "DHITsigNeg",
                                                    DHITsig_PRPS_class == "DHITsigPos" ~ "DHITsigPos",
                                                    DHITsig_PRPS_class == "DHITsig+" ~ "DHITsigPos",
                                                    DHITsig_PRPS_class == "DHITsigNeg" ~ "DHITsigNeg",
                                                    DHITsig_PRPS_class == "DHITsig-" ~ "DHITsigNeg",
                                                    DHITsig_PRPS_class == "UNCLASS" ~ "DHITsigPos",
                                                    TRUE ~ "NA"))

  all_meta = all_meta %>%
    dplyr::mutate(DHITsig_consensus = case_when(DHITsig_consensus == "NEG" ~ "DHITsigNeg",
                                                DHITsig_consensus == "POS" ~ "DHITsigPos",
                                                DHITsig_consensus == "UNCLASS" ~ "DHITsig-IND",
                                                DHITsig_consensus == "DHITsigNeg" ~ DHITsig_consensus,
                                                DHITsig_consensus == "DHITsigPos" ~ DHITsig_consensus,
                                                TRUE ~ "NA"))

  #assign a rank to each pathology for consistent and sensible ordering
  all_meta = all_meta %>%
    dplyr::mutate(pathology_rank = case_when(pathology == "B-ALL" ~ 0,
                                             pathology == "SCBC" ~ 2,
                                             pathology == "CLL" ~ 3,
                                             pathology == "MCL" ~ 4,
                                             pathology == "BL" ~ 7,
                                             pathology == "DLBCL-BL-like" ~ 10,
                                             pathology == "HGBL" ~ 11,
                                             pathology == "COMFL" ~ 13,
                                             pathology == "FL" ~ 15,
                                             pathology == "DLBCL" ~ 19,
                                             pathology == "PBL" ~ 27,
                                             pathology == "B-cell unclassified" ~ 29,
                                             pathology == "MM" ~ 33,
                                             TRUE ~ 35))

  all_meta = all_meta %>%
    dplyr::mutate(lymphgen_rank = case_when(pathology != "DLBCL" ~ pathology_rank,
                                            lymphgen == "Other" ~ 16,
                                            lymphgen == "COMPOSITE" ~ 17,
                                            lymphgen == "N1" ~ 18,
                                            lymphgen == "EZB" ~ 19,
                                            lymphgen == "ST2" ~ 20,
                                            lymphgen == "BN2" ~ 21,
                                            lymphgen == "MCD" ~ 22,
                                            TRUE ~ 50))
  if(with_outcomes){

    all_meta = left_join(all_meta, outcome_table, by = "patient_id") %>%
      mutate(age_group = case_when(cohort == "BL_Adult"~"Adult_BL", cohort == "BL_Pediatric" | cohort == "BL_ICGC" ~ "BL_Pediatric", TRUE ~ "Other"))

  }
  # take one row per sample_id using seq_type_priority
  if(seq_type_priority=="genome"){
    all_meta = all_meta %>%
      arrange(sample_id,seq_type) %>%
      group_by(sample_id) %>%
      slice_tail() %>%
      ungroup()
  }
  if(seq_type_priority=="capture"){
    all_meta = all_meta %>%
      arrange(sample_id,seq_type) %>%
      group_by(sample_id) %>%
      slice_head() %>%
      ungroup()
  }
  return(all_meta)
}
