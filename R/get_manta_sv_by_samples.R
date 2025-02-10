#' @title Get Manta SV By Samples.
#'
#' @description Load the manta output for a set of samples.
#'
#' @details This is a convenience wrapper function for [GAMBLR.results::get_manta_sv_by_sample] (and called by [GAMBLR.results::get_manta_sv]).
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_combined_sv],
#' [GAMBLR.results::get_manta_sv], [GAMBLR.results::get_manta_sv_by_sample]
#'
#' @param these_samples_metadata The only required parameter is a metadata table (data frame) that must contain a row for each sample you want the data from. The additional columns the data frame needs to contain, besides sample_id, are: unix_group, genome_build, seq_type, pairing_status.
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Default value is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned. Default value is 40.
#' @param pass If set to TRUE, only return SVs that are annotated with PASS in the FILTER column. Set to FALSE to keep all variants, regardless if they PASS the filters. Default is TRUE.
#' @param projection The projection of returned calls. Default is grch37.
#' @param verbose Set to FALSE to prevent the path of the requested bedpe file to be printed.
#'
#' @return A data frame containing the Manta outputs from all sample_id in these_samples_metadata in a bedpe-like format with additional columns extracted from the VCF column.
#'
#' @import dplyr stringr
#' @export
#'
#' @examples
#' all_sv = get_manta_sv()
#' metadata = get_gambl_metadata()
#' missing_samples = dplyr::anti_join(metadata, all_sv, by = c("sample_id" = "tumour_sample_id"))
#'
#' missing_from_merge = get_manta_sv_by_samples(these_samples_metadata = missing_samples, verbose = FALSE)
#'
#' @keywords internal
get_manta_sv_by_samples = function(these_samples_metadata,
                                   min_vaf = 0.1,
                                   min_score = 40,
                                   pass = TRUE,
                                   projection = "grch37",
                                   verbose = TRUE){

  #check remote configuration
  remote_session = check_remote_configuration(auto_connect = TRUE)

  #get sample IDs from metadata.
  samples = pull(these_samples_metadata, sample_id)

  #create an empty list.
  all_bedpe = list()

  #wrap get_manta_sv_by_sample.
  all_bedpe = lapply(samples, function(x){get_manta_sv_by_sample(this_sample_id = x,
                                                                 these_samples_metadata = these_samples_metadata,
                                                                 force_lift = FALSE, #the wrapper function performs liftover on all samples that need it.
                                                                 return_anyway = TRUE, #make sure unlifted calls, with the extra column (need_lift) are returned.
                                                                 min_vaf = min_vaf,
                                                                 min_score = min_score,
                                                                 pass = pass,
                                                                 projection = projection,
                                                                 verbose = verbose)})

  #un-nest list into long format.
  merged_bedpe = bind_rows(all_bedpe)

  #take out all calls that need to be lifted
  to_be_lifted = merged_bedpe %>%
    dplyr::filter(need_lift == TRUE)

  if(nrow(to_be_lifted) > 0){
    lifted_calls = liftover(data_df = to_be_lifted, target_build = projection)

    #subset calls that does not need a "lift"
    no_lift_needed = merged_bedpe %>%
      dplyr::filter(need_lift == FALSE)

    #combine calls (lifted and not lifted), arrange and sort accordingly, drop temporary column
    merged_bedpe = bind_rows(lifted_calls, no_lift_needed) %>%
      dplyr::select(-need_lift)

  }else{ #if no samples needs to be lifted, just remove the need_lift column
    merged_bedpe = merged_bedpe %>%
      dplyr::select(-need_lift)
  }

  #add chr prefix to the chromosome name for builds that expect it, but only add when necessary
  #hg38 and hg19 (with chr prefix)
  if(projection %in% c("hg38", "hg19")){
    merged_bedpe = merged_bedpe %>%
      dplyr::mutate(CHROM_A = case_when(str_detect(CHROM_A, "chr") ~ CHROM_A, TRUE ~ paste0("chr", CHROM_A))) %>%
      dplyr::mutate(CHROM_B = case_when(str_detect(CHROM_B, "chr") ~ CHROM_B, TRUE ~ paste0("chr", CHROM_B)))
  }

  #grch37 and grch38 (no chr prefix)
  if(projection %in% c("grch37", "grch38")){
    merged_bedpe = merged_bedpe %>%
      dplyr::mutate(CHROM_A = gsub("chr", "", CHROM_A)) %>%
      dplyr::mutate(CHROM_B = gsub("chr", "", CHROM_B))
  }

  #sort data frame
  merged_bedpe = merged_bedpe %>%
    arrange(CHROM_A, CHROM_B, START_A)

  #return merged manta SVs.
  return(merged_bedpe)
}
