#' @title Get GAMBL metadata.
#'
#' @description Return metadata for a selection of samples.
#'
#' @details This function returns metadata for GAMBL samples. Options for subset and filter the returned data are available.
#' For more information on how to use this function with different filtering criteria, refer to the parameter descriptions,
#' examples and vignettes. Embargoed cases (current options: 'BLGSP-study', 'FL-study', 'DLBCL-study', 'FL-DLBCL-study', 'FL-DLBCL-all', 'DLBCL-unembargoed', 'BL-DLBCL-manuscript', 'MCL','MCL-CLL')
#'
#' @param seq_type_filter Filtering criteria (default: all genomes).
#' @param tissue_status_filter Filtering criteria for tissue status. Possible values are "tumour" (the default) or "normal".
#' @param case_set Optional short name for a pre-defined set of cases avoiding any embargoed cases (current options: 'BLGSP-study', 'FL-study', 'DLBCL-study', 'FL-DLBCL-study', 'FL-DLBCL-all', 'DLBCL-unembargoed', 'BL-DLBCL-manuscript', 'MCL','MCL-CLL').
#' @param remove_benchmarking By default the FFPE benchmarking duplicate samples will be dropped.
#' @param sample_flatfile Optionally provide the full path to a sample table to use instead of the default.
#' @param biopsy_flatfile Optionally provide the full path to a biopsy table to use instead of the default.
#' @param with_outcomes Optionally join to gambl outcome data.
#' @param only_available If TRUE, will remove samples with FALSE or NA in the bam_available column (default: TRUE).
#' @param from_flatfile New default is to use the metadata in the flat-files from your clone of the repo. Can be overridden to use the database.
#' @param seq_type_priority For duplicate sample_id with different seq_type available, the metadata will prioritize this seq_type and drop the others. Possible values are "genome" or "capture".
#'
#' @return A data frame with metadata for each biopsy in GAMBL
#' 
#' \describe{
#'   \item{compression}{Format of the original data used as input for our analysis pipelines (cram, bam or fastq)}
#'   \item{bam_available}{Whether or not this file was available when last checked.}
#'   \item{patient_id}{The anonymized unique identifier for this patient. For BC samples, this will be Res ID.}
#'   \item{sample_id}{A unique identifier for the sample analyzed.}
#'   \item{seq_type}{The assay type used to produce this data (one of "genome","capture, "mrna", "promethION")}
#'   \item{genome_build}{The name of the genome reference the data were aligned to.}
#'   \item{tissue_status}{Whether the sample was atumour or normal.}
#'   \item{cohort}{Name for a group of samples that were added together (usually from a single study), often in the format {pathology}_{cohort_descriptor}.}
#'   \item{library_id}{The unique identifier for the sequencing library.}
#'   \item{pathology}{The diagnosis or pathology for the sample}
#'   \item{time_point}{Timing of biopsy in increasing alphabetical order (A = diagnosis, B = first relapse etc)}
#'   \item{protocol}{General protocol for library construction. e.g. "Ribodepletion", "PolyA", or "Genome"}
#'   \item{ffpe_or_frozen}{Whether the nucleic acids were extracted from a frozen or FFPE sample}
#'   \item{read_length}{The length of reads (required for RNA-seq libraries)}
#'   \item{strandedness}{Whether the RNA-seq librayr construction was strand-specific and, if so, which strand. Required for RNAseq; "positive", "negative", or "unstranded")}
#'   \item{seq_source_type}{Required for RNAseq. Usually the same value as ffpe_or_frozen but sometimes immunotube or sorted cells}
#'   \item{data_path}{Symbolic link to the bam or cram file (not usually relevant for GAMBLR)}
#'   \item{link_name}{Standardized naming for symbolic link (not usually relevant for GAMBLR)}
#'   \item{data_path}{Symbolic link to the fastq file (not usually relevant for GAMBLR)}
#'   \item{fastq_link_name}{Standardized naming for symbolic link for FASTQ file, if used (not usually relevant for GAMBLR)}
#'   \item{unix_group}{Whether a source is external and restricted by data access agreements (icgc_dart) or internal (gambl)}
#'   \item{COO_consensus}{TODO}
#'   \item{DHITsig_consensus}{TODO}
#'   \item{COO_PRPS_class}{TODO}
#'   \item{DHITsig_PRPS_class}{TODO}
#'   \item{DLBCL90_dlbcl_call}{TODO}
#'   \item{DLBCL90_dhitsig_call}{TODO}
#'   \item{res_id}{duplicate of sample_id for local samples and NA otherwise}
#'   \item{DLBCL90_code_set}{Code set used for DLBCL90 call. One of DLBCL90 DLBCL90v2 DLBCL90v3}
#'   \item{DLBCL90_dlbcl_score}{TODO}
#'   \item{DLBCL90_pmbl_score}{TODO}
#'   \item{DLBCL90_pmbl_call}{TODO}
#'   \item{DLBCL90_dhitsig_score}{TODO}
#'   \item{myc_ba}{Result from breakapart FISH for MYC locus}
#'   \item{myc_cn}{Result from copy number FISH for MYC locus}
#'   \item{bcl2_ba}{Result from breakapart FISH for BCL2 locus}
#'   \item{bcl2_cn}{Result from copy number FISH for BCL2 locus}
#'   \item{bcl6_ba}{Result from breakapart FISH for BCL6 locus}
#'   \item{bcl6_cn}{Result from copy number FISH for BCL6 locus}
#'   \item{time_since_diagnosis_years}{TODO}
#'   \item{relapse_timing}{TODO}
#'   \item{dtbx}{TODO. OR REMOVE?}
#'   \item{dtdx}{TODO. OR REMOVE?}
#'   \item{lymphgen_no_cnv}{TODO}
#'   \item{lymphgen_with_cnv}{TODO}
#'   \item{lymphgen_cnv_noA53}{TODO}
#'   \item{lymphgen_wright}{The LymphGen call for this sample from Wright et all (if applicable)}
#'   \item{fl_grade}{TODO}
#'   \item{capture_frozen_sample_id}{TODO}
#'   \item{capture_FFPE_sample_id}{TODO}
#'   \item{capture_unknown_sample_id}{TODO}
#'   \item{genome_frozen_sample_id}{TODO}
#'   \item{genome_ctDNA_sample_id}{TODO}
#'   \item{genome_FFPE_sample_id}{TODO}
#'   \item{mrna_PolyA_frozen_sample_id}{TODO}
#'   \item{mrna_Ribodepletion_frozen_sample_id}{TODO}
#'   \item{mrna_Ribodepletion_frozen_sample_id}{TODO}
#'   \item{XXX_cohort}{Cohort name for batch effect correction(?)}
#'   \item{transformation}{TODO}
#'   \item{relapse}{TODO}
#'   \item{ighv_mutation_original}{TODO}
#'   \item{normal_sample_id}{TODO}
#'   \item{pairing_status}{TODO}
#'   \item{ICGC_ID}{TODO}
#'   \item{ICGC_XXX}{metadata value for ICGC cohort inferred from external metadata}
#'   \item{detailed_pathology}{TODO}
#'   \item{COO_final}{TODO}
#'   \item{consensus_pathology}{TODO}
#'   \item{lymphgen}{TODO}
#'   \item{Tumor_Sample_Barcode}{Duplicate of sample_id for simplifying joins to MAF data frames}
#'   \item{consensus_coo_dhitsig}{TODO}
#'   \item{pathology_rank}{Numeric rank for consistent ordering of samples by pathology}
#'   \item{lymphgen_rank}{Numeric rank for consistent ordering of samples by LymphGen}
#'   \item{hiv_status}{TODO}
#'   \item{CODE_XXX}{Event-free status at last follow-up for overall survival (OS), progression-free survival (PFS) etc. 0 = no event/censored. 1 = event}
#'   \item{XXX_YEARS}{Time, in years, from diagnosis to last follow-up for overall survival (OS), progression-free survival (PFS) }
#'   \item{alive}{Theoretically redundant with CODE_OS}
#'   \item{is_adult}{Adult or pediatric at diagnosis. One of "Adult" for adults and "Pediatric" otherwise}
#'   \item{age_group}{Adult_BL or Pediatric_BL or Other, specific to the BLGSP study}
#'   \item{age}{patient age at diagnosis}
#'   \item{sex}{The biological sex of the patient, if available. Allowable options: M, F, NA}
#'   \item{tx_primary}{TODO}
#'   \item{cause_of_death}{TODO}
#' }
#'
#' @import config dplyr tidyr readr RMariaDB DBI GAMBLR.helpers GAMBLR.utils
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

  # check arguments
  stopifnot('`tissue_status_filter` argument should be either "tumour" or "normal".' = {
    identical(tissue_status_filter, "tumour") | identical(tissue_status_filter, "normal")
  })
  
  stopifnot('`seq_type_priority` argument should be either "genome" or "capture".' = {
    identical(seq_type_priority, "genome") | identical(seq_type_priority, "capture")
  })
  

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
    db = GAMBLR.helpers::check_config_value(config::get("database_name"))
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
    case_set_path = GAMBLR.helpers::check_config_value(config::get("sample_sets")$default)
    full_case_set_path =  paste0(GAMBLR.helpers::check_config_value(config::get("repo_base")), case_set_path)
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

  all_meta = GAMBLR.utils::tidy_lymphgen(all_meta,
              lymphgen_column_in = "lymphgen_cnv_noA53",
              lymphgen_column_out = "lymphgen",
              relevel=TRUE)

  #all_meta = GAMBLR.results::collate_lymphgen(all_meta, verbose=FALSE)

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
  
  
  ### For rows that show the same sample_id, always take one that seq_type is "mrna"
  ### and one between "genome" and "capture", based on seq_type_priority. For each 
  ### seq_type, take the row that shows lower count of NAs.
  
  # count NAs per row
  all_meta <- apply( all_meta, 1, function(x) sum(is.na(x)) ) %>% 
    dplyr::mutate(all_meta, na_count = .)
  
  # take `seq_type_priority` lines
  all_meta_res <- dplyr::filter(all_meta, seq_type == seq_type_priority) %>% 
    dplyr::group_by(sample_id) %>% 
    dplyr::slice(which.min(na_count)) %>% 
    dplyr::ungroup()
  samples_with_priority_seq_type <- all_meta_res$sample_id
  
  # take mrna lines (always)
  all_meta_res <- dplyr::filter(all_meta, seq_type == "mrna") %>% 
    dplyr::group_by(sample_id) %>% 
    dplyr::slice(which.min(na_count)) %>% 
    dplyr::ungroup() %>% 
    rbind(all_meta_res, .)
  
  # define which seq type is not priority
  seq_type_non_priotiry <- ifelse(seq_type_priority == "genome", "capture", "genome")
  # take only those lines that seq_type is not `seq_type_priority` (those were already solved)
  all_meta <- dplyr::filter(all_meta, ! sample_id %in% samples_with_priority_seq_type)
  # take non-`seq_type_priority` lines
  all_meta_res <- dplyr::filter(all_meta, seq_type == seq_type_non_priotiry) %>% 
    dplyr::group_by(sample_id) %>% 
    dplyr::slice(which.min(na_count)) %>% 
    dplyr::ungroup() %>% 
    rbind(all_meta_res) %>% 
    dplyr::select(-na_count) %>% 
    dplyr::arrange(sample_id, seq_type)
  
  all_meta_res
}
