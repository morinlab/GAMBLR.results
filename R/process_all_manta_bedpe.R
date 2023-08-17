#' @title Process All Manta Bedpe.
#'
#' @description This function is in draft mode.
#'
#' @details This is a helper function that is not meant to be used routinely.
#'
#' @param file_df Paths to bedpe.
#' @param out_dir output directory.
#' @param group The unix group.
#' @param genome_build Genome build.
#' @param projection_build The genome we want all results to be relative to (lifted if necessary).
#'
#' @import dplyr readr
#' 
#' @noRd
#'
process_all_manta_bedpe = function(file_df,
                                   out_dir,
                                   group,
                                   genome_build,
                                   projection_build = "grch37"){

  to_merge = list()
  if(missing(out_dir)){
    project_base = check_config_value(config::get("project_base"))
    base_out_dir = check_config_value(config::get("results_staging")$manta)
    out_dir = paste0(project_base, group, "/", base_out_dir)
  }

  process_manta = function(bedpe_file, liftover_to_hg19 = FALSE, liftover_to_hg38 = FALSE, only_return_missing = FALSE, projection = "grch37"){
    cnames = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "NAME", "SOMATIC_SCORE", "STRAND_A", "STRAND_B", "TYPE", "FILTER", "VAF_tumour", "VAF_normal", "DP_tumour", "DP_normal", "tumour_sample_id", "normal_sample_id", "pair_status")
    svbed = suppressMessages(read_tsv(bedpe_file, comment = "##", col_types = "cddcddccccccccccccccccc"))
    this_patient = colnames(svbed)[23]
    this_normal = colnames(svbed)[22]

    if(grepl("--unmatched", bedpe_file)){
      pair_status = "unmatched"
      svbed$pair_status = "unmatched"
    }else{
      svbed$pair_status = "matched"
      pair_status = "matched"
    }
    is_lifted = "native"
    if(liftover_to_hg19 || liftover_to_hg38){
      is_lifted = "lifted"
    }

    if(genome_build == projection | (genome_build == "hs37d5" & projection == "grch37")){
      is_lifted = "native"
    }
    out_file = paste0(out_dir, "/", this_patient, "--", this_normal, "--", pair_status, "--", is_lifted, "--genome--", genome_build, "--", projection_build, "_sv.tsv")
    message("working on OVER HERE:", bedpe_file)
    print(paste("output:", out_file))
    if(file.exists(out_file)){
      if(!only_return_missing){
        print(paste("LOADING", out_file))
        svbed = suppressMessages(read_tsv(out_file, col_types = "ccccccccccccnnnnccc", col_names = cnames))
        return(svbed)
      }
      else{
        svbed = dplyr::filter(svbed, is.na(tumour_sample_id))
        return(svbed)
      }
    }
    if(liftover_to_hg19){
      svbed = liftover_bedpe(bedpe_df = svbed)
    }else if(liftover_to_hg38){
      svbed = liftover_bedpe(bedpe_df = svbed, target_build = "hg38")
    }

    infos = pull(svbed, this_patient)
    infos_n = pull(svbed, this_normal)
    colnames(svbed)[c(1:6)] = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B")

    svbed$VAF_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
    svbed$DP_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 2)[1])})
    svbed$VAF_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
    svbed$DP_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 2)[1])})
    svbed$SOMATIC_SCORE = sapply(svbed$INFO_A, function(x){as.numeric(tail(unlist(strsplit(x, "=")), 1))})

    svbed$tumour_sample_id = this_patient
    svbed$normal_sample_id = this_normal
    message(paste("checking status:", bedpe_file))

    svbed$NAME = "."
    svbed = svbed %>%
      dplyr::select(CHROM_A, START_A, END_A, CHROM_B, START_B, END_B, NAME, SOMATIC_SCORE, STRAND_A, STRAND_B, TYPE, FILTER, VAF_tumour, VAF_normal, DP_tumour, DP_normal, tumour_sample_id, normal_sample_id, pair_status)

    #remove chr prefix from both chromosome names
    svbed = svbed %>%
      mutate(CHROM_A = gsub("chr", "", CHROM_A)) %>%
      mutate(CHROM_B = gsub("chr", "", CHROM_B))
    #print(paste("writing output to",out_file))
    #run liftover after formatting?

    write_tsv(svbed, out_file, col_names = FALSE)
    return(svbed)
  }

  #separately run the hg38 and other builds, separately run per unix_group
  if(projection_build == "grch37"){
    if(genome_build == "hg38"){
      hg38_files = dplyr::filter(file_df, genome_build == "hg38" & unix_group == group) %>%
        pull(file_path)

      bed_data_lifted = hg38_files %>%
        purrr::map(process_manta, liftover_to_hg19 = TRUE, projection = projection_build) %>%
        purrr::reduce(rbind)
    }else{
      not_hg38_files = dplyr::filter(file_df, genome_build != "hg38" & unix_group == group) %>%
        pull(file_path)

      bed_data_not_lifted = not_hg38_files %>%
        purrr::map(process_manta, liftover_to_hg19 = FALSE, projection = projection_build) %>%
        purrr::reduce(rbind)
    }
  }else if(projection_build == "hg38"){
    if(genome_build == "hg38"){
      hg38_files = dplyr::filter(file_df, genome_build == "hg38" & unix_group == group) %>%
        pull(file_path)

      bed_data_not_lifted = hg38_files %>%
        purrr::map(process_manta, liftover_to_hg38 =FALSE, liftover_to_hg19 = FALSE, projection = projection_build) %>%
        purrr::reduce(rbind)
    }else{
      not_hg38_files = dplyr::filter(file_df, genome_build != "hg38" & unix_group == group) %>%
        pull(file_path)

      bed_data_lifted = not_hg38_files %>%
        purrr::map(process_manta, liftover_to_hg38 = TRUE, liftover_to_hg19 = FALSE, projection = projection_build) %>%
        purrr::reduce(rbind)
    }
  }
}
