#' @title Get GAMBL metadata.
#'
#' @description Return metadata for a selection of samples.
#'
#' @details This function returns metadata for GAMBL samples. 
#' This replaces the functionality of the original version, which is still available under the new name [GAMBLR.results::og_get_gambl_metadata]. 
#' The purpose of this function is to provide the metadata for the non-redundant set of samples from GAMBL, dealing with all 
#' types of redundancy caused by samples or biopsies that have data from >1 seq_type (genome or capture), different capture protocols (exome or targeted capture) etc.
#' 
#' @param dna_seq_type_priority The default is "genome" and the only other option is "capture". For duplicate biopsy_id/patient combinations with different seq_type available, prioritize this seq_type and drop the others. 
#' @param capture_protocol_priority For duplicate biopsy_id/patient combinations with different seq_type available, prioritize this seq_type and drop the others. 
#' @param exome_capture_space_priority A vector specifying how to prioritize exome capture space #TODO: implement and test once examples are available
#' @param mrna_collapse_redundancy Default: TRUE. Set to FALSE to obtain all rows for the mrna seq_type including those that would otherwise be collapsed. 
#' @param also_normals Set to TRUE to force the return of rows where tissue_status is normal (default is to restrict to tumour) 
#' @param invert Set to TRUE to force the function to return only the rows that are lost in all the prioritization steps (mostly for debugging)
#' @param everything Set to TRUE to include samples with `bam_available == FALSE`. Default: FALSE - only samples with `bam_available = TRUE` are retained.
#' @param verbose Set to TRUE for a chatty output (mostly for debugging)
#' 
#' @return A data frame with metadata for each biopsy in GAMBL
#' 
#' \describe{
#'   \item{compression}{Format of the original data used as input for our analysis pipelines (cram, bam or fastq)}
#'   \item{bam_available}{Whether or not this file was available when last checked.}
#'   \item{patient_id}{The anonymized unique identifier for this patient. For BC samples, this will be Res ID.}
#'   \item{sample_id}{A unique identifier for the sample analyzed.}
#'   \item{seq_type}{The assay type used to produce this data (one of "genome","capture, "mrna", "promethION")}
#'   \item{capture_space}{Unique ID for the capture space, where applicable}
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
#'   \item{fastq_data_path}{Symbolic link to the fastq file (not usually relevant for GAMBLR)}
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
#'   \item{transformation}{TODO}
#'   \item{relapse}{TODO}
#'   \item{ighv_mutation_original}{TODO}
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
#' tumour_and_normal_metadata = get_gambl_metadata(tissue_status_filter = c('tumour','normal'))
#'
#' non_redundant_genome_and_capture = get_gambl_metadata(seq_type_filter = c('genome', 'capture'),
#'                                                        seq_type_priority = "genome")
#'                                                        
#' absolutely_everything = get_gambl_metadata(seq_type_filter = c('genome', 'capture','mrna'), tissue_status_filter=c('tumour','normal'))
#'
get_gambl_metadata = function(dna_seq_type_priority = "genome", 
                               capture_protocol_priority = "Exome", 
                               dna_preservation_priority = "frozen",
                               exome_capture_space_priority = c("agilent-sureselect-human-all-exon-v7",
                                                                "agilent-sureselect-V5-plus-utr",
                                                                "idt-xgen-v2-grch37",
                                                                "exome-utr-grch37",
                                                                "exome-utr-grch38",
                                                                "none"),
                               mrna_collapse_redundancy=TRUE, 
                               also_normals = FALSE,
                               everything=FALSE,
                               verbose=FALSE,
                               invert=FALSE,
                              ...){
  if(any(names(match.call(expand.dots = TRUE)) %in% formalArgs(og_get_gambl_metadata))){
    args_match = names(match.call(expand.dots = TRUE))[which(names(match.call(expand.dots = TRUE)) %in% formalArgs(og_get_gambl_metadata))]
    message("one or more arguments for the original get_gambl_metadata detected, reverting to that function")
    print(args_match)
    return(og_get_gambl_metadata(...))
  }
  dropped_rows = list()
  check_remote_configuration()
  #this needs to be in any function that reads files from the bundled GAMBL outputs synced by Snakemake
  
  base = config::get("repo_base")
  sample_flatfile = paste0(base, config::get("table_flatfiles")$samples)
  sample_meta = suppressMessages(read_tsv(sample_flatfile, guess_max = 100000))
  
  
  
  biopsy_flatfile = paste0(base, config::get("table_flatfiles")$biopsies)
  biopsy_meta = suppressMessages(read_tsv(biopsy_flatfile, guess_max = 100000))
  
  
  
  # We virtually always want to remove samples without bam_available == TRUE
  if(everything){
    mrna_collapse_redundancy = FALSE
  }else{
    sample_meta = dplyr::filter(sample_meta, bam_available %in% c(1, "TRUE"))
  }
  
  sample_meta_normal =  sample_meta %>%
    dplyr::filter(tissue_status == "normal") 

  sample_meta_tumour =  sample_meta %>%
    dplyr::filter(tissue_status != "normal") 
  
  massage_tumour_metadata = function(tumour_metadata){
    #check that capture samples have a protocol, fill in missing values and warn about it
    num_missing_protocol = filter(tumour_metadata,seq_type=="capture",is.na(protocol)) %>% nrow()
    message(paste(num_missing_protocol,"capture samples are missing a value for protocol. Assuming Exome."))
    tumour_metadata = mutate(tumour_metadata,protocol=case_when(seq_type == "capture" & is.na(protocol) ~ "Exome",
                                                                seq_type == "genome" & is.na(protocol) ~ "Genome",
                                                                TRUE ~ protocol))
    return(tumour_metadata)
  }
  sample_meta_tumour = massage_tumour_metadata(sample_meta_tumour)
  
  #currently this function just nags the user
  check_biopsy_metadata = function(tumour_metadata){
    base = config::get("repo_base")
    flatfile = paste0(base, config::get("table_flatfiles")$biopsies)
    b_meta = suppressMessages(read_tsv(flatfile, guess_max = 100000))
    #sanity check biopsy_metadata contents
    missing_biopsies = filter(tumour_metadata,!biopsy_id %in% b_meta$biopsy_id) %>% select(sample_id,biopsy_id,cohort,pathology)
    n_missing_biopsies = nrow(missing_biopsies)
    if(n_missing_biopsies>0){
      message(paste(n_missing_biopsies,"biopsies are missing from the biopsy metadata. This should be fixed!"))
      cohorts_missing = unique(missing_biopsies$cohort)
      message(paste("affected cohorts: ", paste(cohorts_missing,collapse=",")))
    }
    #compare the shared columns for consistency and warn about any inconsistencies
    combined_meta = left_join(tumour_metadata,b_meta,by=c("patient_id","biopsy_id"))
    pathology_discrepancies = filter(combined_meta,pathology.x != pathology.y)
    num_path = nrow(pathology_discrepancies)
    time_point_discrepancies = filter(combined_meta,time_point.x != time_point.y)
    num_time_point = nrow(time_point_discrepancies)
    if(num_path>0){
      message(paste(num_path, "biopsies with discrepancies in the pathology field. This should be fixed!"))
      if(verbose){
        dplyr::select(pathology_discrepancies,biopsy_id,patient_id,pathology.x,pathology.y) %>% print()
      }
    }
    if(num_time_point>0){
      message(paste(num_time_point, "biopsies with discrepancies in the time_point field. This should be fixed!"))
      if(verbose){
        dplyr::select(time_point_discrepancies,biopsy_id,patient_id,time_point.x,time_point.y) %>% print()
      }
    }
    
  }
  
  check_biopsy_metadata(sample_meta_tumour)
  
  sample_meta_rna = filter(sample_meta,seq_type == "mrna") #includes some normals
  if(mrna_collapse_redundancy){
    sample_subset_rna = check_gene_expression(verbose=verbose) %>% select(-protocol,-cohort,-ffpe_or_frozen)
    sample_meta_rna_dropped = suppressMessages(anti_join(sample_meta_rna,sample_subset_rna))
    sample_meta_rna_kept = suppressMessages(inner_join(sample_meta_rna,sample_subset_rna))
    dropped_rows[["mrna"]] = sample_meta_rna_dropped
  }else{
    sample_meta_rna_kept = sample_meta_rna
  }
  
  sample_meta_tumour_dna = sample_meta_tumour %>%
    dplyr::filter(seq_type %in% c("genome","capture")) 
  
  sample_meta_normal_dna = sample_meta_normal %>% 
    dplyr::filter(seq_type %in% c("genome","capture")) 

  
  #helper function to prioritize systematically on the values in a column specified by column_name
  prioritize_metadata_column = function(metadata,grouping_column_name,priority_column_name,priority_value){
    #ALWAYS group by patient_id, biopsy_id
    original = nrow(metadata)
    if(missing(grouping_column_name)){
      dupes = group_by(metadata, patient_id, biopsy_id) %>% tally() %>% filter(n>1)
    }else{
      dupes = group_by(metadata, patient_id, biopsy_id,{{grouping_column_name}}) %>% tally() %>% filter(n>1)
    }

    if(missing(grouping_column_name)){
      collapsed = mutate(metadata,priority=ifelse({{priority_column_name}}==priority_value,1,2)) %>%
        group_by( patient_id, biopsy_id) %>%
        arrange(priority,sample_id,.by_group = T) %>%  #arrange on sample_id for the rare case of >1 sample ID from the same biopsy
        slice_head(n=1)
    }else{
      collapsed = mutate(metadata,priority=ifelse({{priority_column_name}}==priority_value,1,2)) %>%
        group_by( patient_id, biopsy_id,{{grouping_column_name}}) %>%
        arrange(priority,sample_id,.by_group=T) %>% 
        slice_head(n=1) %>%
        ungroup()
    }
    now = nrow(collapsed)
    if(verbose){
      message(paste("started with ",original,"rows and now have ",now))
    }
    #collapsed = select(collapsed,-priority)
    return(ungroup(collapsed))
  }
  sample_meta_tumour_dna_capture = filter(sample_meta_tumour_dna,seq_type=="capture")
  
  sample_meta_tumour_dna_genome = filter(sample_meta_tumour_dna,seq_type=="genome")
  
  #prioritize within the genome seq_type (drop any FFPE where we have a frozen)
  if(verbose){
    message(paste("prioritizing",dna_preservation_priority))
  }
  sample_meta_tumour_dna_genome_kept = prioritize_metadata_column(sample_meta_tumour_dna_genome,
                                                                  priority_column_name = ffpe_or_frozen,
                                                                  priority_value = dna_preservation_priority) 
  
  sample_meta_tumour_dna_genome_dropped = filter(sample_meta_tumour_dna_genome,!sample_id %in% sample_meta_tumour_dna_genome_kept$sample_id)
  dropped_rows[["tumour_genome"]] = sample_meta_tumour_dna_genome_dropped

  #prioritize within the capture seq_type
  if(verbose){
    message(paste("prioritizing",capture_protocol_priority, "for the capture protocol"))
  }
  sample_meta_tumour_dna_capture_kept = prioritize_metadata_column(sample_meta_tumour_dna_capture,
                                                                   seq_type,
                                                                   protocol,
                                                                   capture_protocol_priority) 

  sample_meta_tumour_dna_capture_dropped = suppressMessages(anti_join(sample_meta_tumour_dna_capture,sample_meta_tumour_dna_capture_kept)) %>% unique()
  dropped_rows[["tumour_capture"]] = sample_meta_tumour_dna_capture_dropped
  
  #prioritize across the DNA seq_types using dna_seq_type_priority paraameter
  sample_meta_tumour_dna = bind_rows(sample_meta_tumour_dna_genome,sample_meta_tumour_dna_capture_kept)
  
  if(verbose){
    message(paste("prioritizing",dna_seq_type_priority, "for tumours"))
  }
  sample_meta_tumour_dna_kept = prioritize_metadata_column(sample_meta_tumour_dna,
                                                           priority_column_name=seq_type,
                                                           priority_value=dna_seq_type_priority) 
  
  sample_meta_tumour_dna_dropped = suppressMessages(anti_join(sample_meta_tumour_dna,sample_meta_tumour_dna_kept))
  dropped_rows[["tumour_dna"]] = sample_meta_tumour_dna_dropped
  
  if(verbose){
    message(paste("prioritizing",dna_seq_type_priority, "for normals"))
  }
  sample_meta_normal_dna_kept = prioritize_metadata_column(sample_meta_normal_dna,
                                                           priority_column_name=seq_type,
                                                           priority_value=dna_seq_type_priority)
  
  sample_meta_normal_dna_dropped = suppressMessages(anti_join(sample_meta_normal_dna,sample_meta_normal_dna_kept))
  dropped_rows[["normal_dna"]] = sample_meta_normal_dna_dropped
  
  
  if(invert){
    #give the user back only the metadata rows that would otherwise be suppressed/dropped by the function
    all_dropped = do.call("bind_rows",dropped_rows) %>% unique()
    if(verbose){
      message("Overview of samples dropped due to redundancy/prioritization")
      group_by(all_dropped,seq_type) %>% tally() %>% print()
    }
    return(all_dropped)
  }
  #remove redundant columns before joining, preferring the values in gambl_samples_available
  biopsy_meta = select(biopsy_meta,-time_point,-pathology,-EBV_status_inf)
  
  
  
  all_meta_kept = bind_rows(sample_meta_tumour_dna_kept, filter(sample_meta_rna_kept,tissue_status=="tumour")) %>% 
    ungroup() %>% select(-priority)
  all_meta_kept = left_join(all_meta_kept,biopsy_meta,by=c("patient_id","biopsy_id")) 
  all_meta_kept = GAMBLR.utils::tidy_lymphgen(all_meta_kept,
                                         lymphgen_column_in = "lymphgen_cnv_noA53",
                                         lymphgen_column_out = "lymphgen",
                                         relevel=TRUE)
  if(also_normals){
    #add missing normals to the data frame by calling the function again
    all_meta_kept = bind_rows(all_meta_kept,sample_meta_normal_dna_kept,filter(sample_meta_rna_kept,tissue_status=="normal"))
    return(all_meta_kept)
  }else{
    return(all_meta_kept)  
  }
  
}
