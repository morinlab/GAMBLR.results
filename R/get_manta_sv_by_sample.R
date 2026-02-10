#' @title Get Manta SV By Sample.
#'
#' @description Load the manta output (from individual flat file) for 1 sample.
#'
#' @details This function is used for retrieving Manta results (structural variants) from individual flat-files (one sample).
#' For multiple samples, please see [GAMBLR.results::get_manta_sv_by_samples] (a convenience wrapper function for [GAMBLR.results::get_manta_sv_by_sample]).
#' Additional columns are extracted from the VCF column and standard filtering options are available.
#' This function also performs a lift-over to selected projection, if needed.
#' Please note, if `force_lift` is set to FALSE, an extra column will be added that states if the returned variant calls need to be lifted.
#' The value for this column is returned TRUE (for all rows) if the available genome projection for the selected sample does not match the selected projection (i.e requiring the user to manually lift the calls).
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_combined_sv], [GAMBLR.results::get_manta_sv], [GAMBLR.results::get_manta_sv_by_samples]
#'
#' @param this_sample_id The single sample ID you want to obtain the result from. If this parameter is not supplied, the function will retrieve sample ID from the supplied metadata table (these_samples_metadata).
#' @param these_samples_metadata A metadata table containing metadata for this_sample_id, or sample of interest. This parameter is required.
#' @param force_lift If TRUE, coordinates will be lifted (if needed) to the selected projection. Default is FALSE. WARNING: if your code calls this function directly, set this parameter to TRUE to ensure that the returned calls are in respect to the requested projection.
#' @param return_anyway Set to TRUE to force variant calls to be returned, even if they're not lifted, This parameter should only ever be modified from the default setting when this function is called by another function that handles the liftOver separately.
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Default value is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned. Default value is 40.
#' @param pass_filters If set to TRUE, only return SVs that are annotated with PASS in the FILTER column. Set to FALSE to keep all variants, regardless if they PASS the filters. Default is TRUE.
#' @param projection The projection of returned calls. Default is grch37.
#' @param verbose Set to FALSE to prevent the path of the requested bedpe file to be printed.
#'
#' @return a data frame containing the Manta outputs from this_sample_id in a bedpe-like format with additional columns extracted from the VCF column.
#'
#' @import config dplyr readr stringr tibble glue GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' #example 1
#' #get manta calls for a sample that needs to be lifted to "hg38" and let this function
#' #take care of the liftover step for you.
#' my_sv = get_manta_sv_by_sample(this_sample_id = "99-27783_tumorA",
#'                                these_samples_metadata = get_gambl_metadata(),
#'                                projection = "hg38",
#'                                force_lift = TRUE)
#'
#' #example 2
#' #get manta calls based on an already filtered metadata (with one sample ID)
#' my_metadata = get_gambl_metadata()
#' my_metadata = dplyr::filter(my_metadata, sample_id=="99-27783_tumorA")
#'
#' my_sv = get_manta_sv_by_sample(these_samples_metadata = my_metadata,
#'                                projection = "hg38",
#'                                force_lift = TRUE)
#' }
#' @keywords internal
get_manta_sv_by_sample = function(this_sample_id,
                                  these_samples_metadata,
                                  force_lift = FALSE,
                                  return_anyway = FALSE,
                                  min_vaf = 0.1,
                                  min_score = 40,
                                  pass_filters = TRUE,
                                  projection = "grch37",
                                  verbose = TRUE){

  #safetynet for preventing users to mistakenly return un-lifted variant calls.
  if(!force_lift){ #i.e I will run liftover on my own, based on the information in the extra column (need_lift).
    if(!return_anyway){
      stop("If you know what you are doing and wish to liftover the returned sample yourself, set return_anyway to TRUE. If you want this function to handle the liftover for you, set force_lift = TRUE")
    }
  }

  #check remote configuration
  remote_session = check_remote_configuration(auto_connect = TRUE)

  if(missing(this_sample_id)){
    if(!nrow(these_samples_metadata) == 1){
      stop("There is more than one sample in the supplied metadata table. Either subset metadata to only have one sample, provide the this_sample_id parameter OR consider running get_manta_sv_by_samples")
    }
    this_sample_id = these_samples_metadata$sample_id
  }

  these_samples_metadata = dplyr::filter(these_samples_metadata, sample_id == this_sample_id)
  if(!nrow(these_samples_metadata==1)){
    stop("metadata does not seem to contain your this_sample_id or you didn't provide one")
  }

  #get wildcards
  tumour_sample_id = this_sample_id
  unix_group = pull(these_samples_metadata, unix_group)
  seq_type = pull(these_samples_metadata, seq_type)
  genome_build = pull(these_samples_metadata, genome_build)
  pairing_status = pull(these_samples_metadata, pairing_status)

  if(pairing_status == "matched"){
    normal_sample_id = pull(these_samples_metadata, normal_sample_id)
  }else{
    normal_sample_id = config::get("unmatched_normal_ids")[[unix_group]][[seq_type]][[genome_build]]
  }

  #get samples from individual flat files
  path_template = GAMBLR.helpers::check_config_value(config::get("results_flatfiles")$sv_manta$template)

  if(!remote_session){
    path_template_full = paste0(GAMBLR.helpers::check_config_value(config::get("project_base")), path_template)
    bedpe_path = glue::glue(path_template_full)
    if(!file.exists(bedpe_path)){
      print(paste("missing: ", bedpe_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }
  }else{
    local_path_template = paste0(GAMBLR.helpers::check_config_value(config::get("project_base", config = "remote")), path_template)
    bedpe_path = glue::glue(local_path_template)

    #check if the requested file is on your local machine, if not, get it!
    if(!file.exists(bedpe_path)){
      remote_path_template = paste0(GAMBLR.helpers::check_config_value(config::get("project_base", config = "default")), path_template)
      remote_bedpe_path = glue::glue(remote_path_template)
      cat(paste0("Local file not found.\ntrying to copy requested file: ", remote_bedpe_path, "\n", "To: ", bedpe_path))
      dirN = dirname(bedpe_path)
      suppressMessages(suppressWarnings(dir.create(dirN, recursive = T)))
      ssh::scp_download(ssh_session, remote_bedpe_path, dirN)
    }
  }

  #read sample flat-file
  if(verbose){
    message(paste0("Reading ", this_sample_id, " from: ", bedpe_path))
  }
  if(file.exists(bedpe_path)){
    bedpe_dat_raw = suppressMessages(read_tsv(bedpe_path,
                                                progress = FALSE,
                                                comment = "##",
                                                col_types = "cddcddccccccccccccccccc"))
  }else{
    message("DID NOT FIND THE FILE AT THE EXPECTED PATH!!!")
    message(paste0("Skipping ", this_sample_id, " from: ", bedpe_path))
    return()
  }

  #return empty data frame
  if(!nrow(bedpe_dat_raw==0)){
    message(paste0("WARNING! No SV calls found in flat-file for: ", this_sample_id))
    return()
  }

  #if the selected projection is different from the genome build (for the selected sample), add information that this sample needs to be lifted (by get_manta_sv_by_samples).
  if(genome_build != projection){
    if(all(str_detect(genome_build, "37|19"))){
      if(projection %in% c("grch37", "hg19")){
        bedpe_dat_raw = bedpe_dat_raw %>%
          add_column(need_lift = FALSE)
      }else if(projection %in% c("hg38", "grch38")){
        bedpe_dat_raw = bedpe_dat_raw %>%
          add_column(need_lift = TRUE)
      }else{
        stop(paste0(projection, " is not a valid projection, acceptable projections are; grch37, grch38, hg19, hg38"))
      }
    }else if(all(str_detect(genome_build, "38"))){
      if(projection %in% c("grch37", "hg19")){
        bedpe_dat_raw = bedpe_dat_raw %>%
          add_column(need_lift = TRUE)
      }else if(projection %in% c("hg38", "grch38")){
        bedpe_dat_raw = bedpe_dat_raw %>%
          add_column(need_lift = FALSE)
      }else{
        stop(paste0(projection, " is not a valid projection, acceptable projections are; grch37, grch38, hg19, hg38"))
      }
    }
  }else{
    bedpe_dat_raw = bedpe_dat_raw %>%
      add_column(need_lift = FALSE)
  }

  if(force_lift){
    if(bedpe_dat_raw$need_lift[1] == TRUE){
      bedpe_dat_raw = liftover(data_df = bedpe_dat_raw, target_build = projection)
      message(paste0(this_sample_id, " flat-file is not available in the selected projection, running liftover_bedpe..."))
      message(paste0(this_sample_id, " successfully lifted to ", projection))
    }
  }

  #data wrangling
  #get infos
  infos = pull(bedpe_dat_raw, tumour_sample_id)
  infos_n = pull(bedpe_dat_raw, normal_sample_id)

  #create new columns with sample IDs
  bedpe_dat = bedpe_dat_raw %>%
    mutate(tumour_sample_id = tumour_sample_id, normal_sample_id = normal_sample_id, pair_status = pairing_status)

  #rename columns to match the expected format
  colnames(bedpe_dat)[c(1:6)] = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B")

  #extract info fields from VCF
  #tumour sample
  bedpe_dat$VAF_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
  bedpe_dat$DP_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")),2)[1])})

  #normal sample
  bedpe_dat$VAF_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
  bedpe_dat$DP_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 2)[1])})

  #get somatic score
  bedpe_dat$SCORE = sapply(bedpe_dat$INFO_A, function(x){as.numeric(tail(unlist(strsplit(x, "=")), 1))})

  #Rename and select columns to match what is returned with get_combined_sv.
  bedpe_dat = bedpe_dat %>%
    rename("DP" = "DP_tumour", "manta_name" = "ID") %>%
    dplyr::select("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B",
                  "manta_name", "SCORE", "STRAND_A", "STRAND_B", "tumour_sample_id",
                  "normal_sample_id", "VAF_tumour", "DP", "pair_status", "FILTER", "need_lift")

  #VAF and somatic score filtering.
  bedpe_dat = bedpe_dat %>%
    dplyr::filter(VAF_tumour >= min_vaf & SCORE >= min_score)

  #Filter on FILTER (variant callers variant filter criteria).
  if(pass_filters){
    bedpe_dat = bedpe_dat %>%
      dplyr::filter(FILTER == "PASS")
  }

  #Deal with chr prefixes based on projection
  if(force_lift){
    #hg38 and hg19 (with chr prefix)
    if(projection %in% c("hg38", "hg19")){
      bedpe_dat = bedpe_dat %>%
        dplyr::mutate(CHROM_A = case_when(str_detect(CHROM_A, "chr") ~ CHROM_A, TRUE ~ paste0("chr", CHROM_A))) %>%
        dplyr::mutate(CHROM_B = case_when(str_detect(CHROM_B, "chr") ~ CHROM_B, TRUE ~ paste0("chr", CHROM_B)))
    }

    #grch37 and grch38 (no chr prefix)
    if(projection %in% c("grch37", "grch38")){
      bedpe_dat = bedpe_dat %>%
        dplyr::mutate(CHROM_A = gsub("chr", "", CHROM_A)) %>%
        dplyr::mutate(CHROM_B = gsub("chr", "", CHROM_B))
    }

    #remove the additional column (need_lift)
    bedpe_dat = bedpe_dat %>%
      dplyr::select(-need_lift)
  }

  #enforce column types and sort returned calls
  bedpe_dat = bedpe_dat %>%
    mutate(across(c(CHROM_A, CHROM_B, manta_name, STRAND_A, STRAND_B, tumour_sample_id, normal_sample_id, pair_status, FILTER), as.character)) %>%
    mutate(across(c(START_A, END_A, START_B, END_B, SCORE, VAF_tumour, DP), as.numeric)) %>%
    arrange(CHROM_A, CHROM_B, START_A)
  bedpe_dat = create_genomic_data(bedpe_dat,projection)
  return(bedpe_dat)
}
