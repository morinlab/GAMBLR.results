#' @title Get SSM By Sample.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample.
#'
#' @details This was implemented to allow flexibility because there are some samples that we may want to use a different set of variants than those in the main GAMBL merge.
#' The current use case is to allow a force_unmatched output to be used to replace the SSMs from the merge for samples with known contamination in the normal.
#' This will also be useful to apply a blacklist to individual MAFs when coupled with [GAMBLR.results::annotate_ssm_blacklist].
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_coding_ssm], [GAMBLR.results::get_coding_ssm_status],
#' [GAMBLR.results::get_ssm_by_patients], [GAMBLR.results::get_ssm_by_samples], [GAMBLR.results::get_ssm_by_region], [GAMBLR.results::get_ssm_by_regions]
#'
#' @param these_samples_metadata Required if not specifying both this_sample_id and this_seq_type a single row or entire metadata table containing your sample_id.
#' @param tool_name The name of the variant calling pipeline (currently only slms-3 is supported).
#' @param projection The projection genome build. Supports hg38 and grch37.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF.
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatibility with additional variant calling parameters/versions.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters.
#' @param verbose Enable for debugging/noisier output.
#' @param this_sample_id Deprecated. Inferred from these_samples_metadata
#' @param this_seq_type Deprecated. Inferred from these_samples_metadata
#'
#' @return data frame in MAF format.
#'
#' @import dplyr tidyr glue GAMBLR.helpers
#'
#' @examples
#' 
#' \dontrun{
#'  this_sample_df = GAMBLR.results:::get_ssm_by_sample(
#'                           these_samples_metadata = get_gambl_metadata() %>%
#'                                 dplyr::filter(sample_id == "HTMCP-01-06-00485-01A-01D",
#'                                      seq_type == "genome"),
#'                          projection = "grch37")
#'
#'  capture_meta = get_gambl_metadata() %>% dplyr::filter(seq_type == "capture")
#'
#'  ssm_sample = GAMBLR.results:::get_ssm_by_sample(this_sample_id = "CASA0002_2015-03-10",
#'                                projection = "grch37",
#'                                augmented = T,
#'                                these_samples_metadata = capture_meta)
#' }
get_ssm_by_sample = function(these_samples_metadata,
                             tool_name = "slms-3",
                             projection = "grch37",
                             augmented = TRUE,
                             flavour = "clustered",
                             min_read_support = 3,
                             basic_columns = TRUE,
                             maf_cols = NULL,
                             verbose = FALSE,
                             this_sample_id,
                             this_seq_type
                             ){
  remote_session = check_remote_configuration(auto_connect = TRUE)
  if(missing(this_sample_id) & missing(these_samples_metadata)){
    stop("Must provide a sample_id or a single row of metadata")
  }else if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata() %>%
      dplyr::filter(sample_id == this_sample_id) %>% 
      dplyr::filter(seq_type != "mrna") 
  }else{
    #just in case more than our row was provided
    if(nrow(these_samples_metadata)>1){
      these_samples_metadata = these_samples_metadata %>%
        dplyr::filter(seq_type != "mrna") %>%
        dplyr::filter(sample_id == this_sample_id)
    }
  }
  if(nrow(these_samples_metadata) > 1 | nrow(these_samples_metadata)==0){
    print(nrow(these_samples_metadata))
    print(these_samples_metadata)
    stop("problem!")
  }
  
  seq_type = pull(these_samples_metadata,seq_type)

  sample_id = pull(these_samples_metadata,sample_id)
  tumour_sample_id = sample_id
  unix_group = pull(these_samples_metadata, unix_group)
  genome_build = pull(these_samples_metadata, genome_build)
  target_builds = projection
  pair_status = pull(these_samples_metadata, pairing_status)
  if(verbose){
    print(paste("group:", unix_group, "genome:", genome_build))
  }
  # Get unmatched normal if necessary. This is done using the unmatched normals that were added to the GAMBLR config.
  # That will need to be kept up to date if/when any new references are added.
  if(pair_status == "unmatched"){
    keys = paste0("unmatched_normal_ids$",unix_group,"$",seq_type,"$",genome_build)
    normal_sample_id = check_config_and_value(keys)
  }else{
    normal_sample_id = pull(these_samples_metadata, normal_sample_id)
  }
  base_path = ""
  if(flavour == "legacy"){
    warning("Access to the old variant calls is not currently supported in this function")
    warning("Use get_ssm_by_samples to access the legacy flavour")
    # To be fixed maybe if we decide it's needed.
    # Implementation of backwards compatability will be a lot harder because of the old naming scheme.
    return()
  }else if(flavour == "clustered"){
    vcf_base_name = "slms-3.final"
    path_template = check_config_and_value(
      "results_flatfiles$ssm$template$clustered$deblacklisted",
      config_name="default"
      )
    path_complete = unname(unlist(glue::glue(path_template)))
    full_maf_path = paste0(
      check_config_and_value("project_base",
      config_name="default"),
      path_complete)
    local_full_maf_path = paste0(
      check_config_and_value("project_base"),
      path_complete)
    if(augmented){
      path_template = check_config_and_value(
        "results_flatfiles$ssm$template$clustered$augmented",
        config_name="default")
      path_complete = unname(unlist(glue::glue(path_template)))
      aug_maf_path = paste0(check_config_and_value("project_base",config_name="default"), path_complete)
      local_aug_maf_path = paste0(check_config_and_value("project_base"), path_complete)
    }
  }else{
    warning("Currently the only flavour available to this function is 'clustered'")
  }
  if(remote_session){
    #check if file exists
    status = ssh::ssh_exec_internal(ssh_session,command=paste("stat",aug_maf_path),error=F)$status

    # first check if we already have a local copy
    # Load data from local copy or get a local copy from the remote path first
    if(status==0){
      if(verbose){
        print(paste("found:",aug_maf_path))
        print(paste("local home:",local_aug_maf_path))
      }
      dirN = dirname(local_aug_maf_path)

      suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
      if(!file.exists(local_aug_maf_path)){
        print(paste("DOWNLOADING:",aug_maf_path))
        ssh::scp_download(ssh_session,aug_maf_path,dirN)
      }

     #check for missingness
     if(!file.exists(local_aug_maf_path)){
      print(paste("missing: ", local_aug_maf_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
     }

      sample_ssm = fread_maf(local_aug_maf_path) %>%
      dplyr::filter(t_alt_count >= min_read_support)
    }else{
      if(verbose){
        print(paste("will use",full_maf_path))
        print(paste("local home:",local_full_maf_path))
      }
      dirN = dirname(local_full_maf_path)

      suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
      if(!file.exists(local_full_maf_path)){
        print(paste("DOWNLOADING:",aug_maf_path))
        ssh::scp_download(ssh_session,full_maf_path,dirN)
      }
      #check for missingness
      if(!file.exists(local_full_maf_path)){
        print(paste("missing: ", local_full_maf_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote"')
      }

      sample_ssm = fread_maf(local_full_maf_path)
    }
  }else if(augmented && file.exists(aug_maf_path)){
    full_maf_path = aug_maf_path

    #check for missingness
    if(!file.exists(full_maf_path)){
      print(paste("missing: ",full_maf_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    sample_ssm = fread_maf(full_maf_path)
    if(min_read_support){
      # drop poorly supported reads but only from augmented MAF
      sample_ssm = dplyr::filter(sample_ssm, t_alt_count >= min_read_support)
    }
  }else{
    if(!file.exists(full_maf_path)){
      print(paste("missing: ", full_maf_path))
      message(paste("warning: file does not exist, skipping it.", full_maf_path))
      return()
    }
    #check for missingness
    if(!file.exists(full_maf_path)){
      print(paste("missing: ", full_maf_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    sample_ssm = fread_maf(full_maf_path)
  }


  #subset maf to only include first 43 columns (default)
  if(basic_columns){
    sample_ssm = dplyr::select(sample_ssm, c(1:45))
    }

  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    sample_ssm = dplyr::select(sample_ssm, all_of(maf_cols))
    }
  sample_ssm = mutate(sample_ssm,maf_seq_type = seq_type)
  sample_ssm = create_maf_data(sample_ssm,projection)
  return(sample_ssm)
}
