#' @title Get SSM By Sample.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample.
#'
#' @details This was implemented to allow flexibility because there are some
#' samples that we may want to use a different set of variants than those in
#' the main GAMBL merge.
#' The current use case is to allow a force_unmatched output to be used
#' to replace the SSMs from the merge for samples with known contamination
#' in the normal.
#' This will also be useful to apply a blacklist to individual MAFs when coupled
#' with [GAMBLR.results::annotate_ssm_blacklist].
#' Is this function not what you are looking for? Try one of the related
#' functions:
#' [GAMBLR.results::get_coding_ssm], [GAMBLR.results::get_coding_ssm_status],
#' [GAMBLR.results::get_ssm_by_samples],
#' [GAMBLR.results::get_ssm_by_region], [GAMBLR.results::get_ssm_by_regions]
#'
#' @param these_samples_metadata Required. A single row of metadata specifying
#' which sample_id and seq_type you desire the mutations from
#' @param this_sample_id Optional. Can be used to subset a multi-row these_samples_metadata 
#' table to a single row. 
#' @param tool_name The name of the variant calling pipeline (currently
#' only slms-3 is supported).
#' @param projection The projection genome build. Supports hg38 and grch37
#' @param augmented default: TRUE. Set to FALSE if you instead want the original
#' MAF from each sample for multi-sample patients instead of the augmented MAF.
#' @param flavour Currently this function only supports one flavour option but
#' this feature is meant for eventual compatibility with additional variant
#' calling parameters/versions.
#' @param min_read_support Only returns variants with at least this many
#' reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 43 columns of MAF rather than full
#' details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify
#' which columns to be returned within the MAF. This parameter can either
#' be a vector of indexes (integer) or a vector of characters.
#' @param verbose Enable for debugging/noisier output.
#' @param this_seq_type Deprecated. Inferred from these_samples_metadata
#'
#' @return data frame in MAF format.
#'
#' @import dplyr tidyr glue GAMBLR.helpers
#'
#' @examples
#' 
#' maf_samp = GAMBLR.results:::get_ssm_by_sample(
#'   get_gambl_metadata() %>% dplyr::filter(sample_id=="13-27975_tumorA"),
#'   augmented = FALSE
#' )
#' nrow(maf_samp)
#' maf_samp_aug = GAMBLR.results:::get_ssm_by_sample(
#'   get_gambl_metadata() %>% dplyr::filter(sample_id=="13-27975_tumorA"),
#'   augmented = TRUE
#' )
#' nrow(maf_samp_aug)
#' 
#' 
#'  some_maf = GAMBLR.results:::get_ssm_by_sample(
#'                           these_samples_metadata = get_gambl_metadata() %>%
#'                             dplyr::filter(sample_id == "HTMCP-01-06-00485-01A-01D",
#'                                      seq_type == "genome"),
#'                          projection = "hg38")
#'  dplyr::select(some_maf,1:10)
#'  
#'
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
  remote_session = check_host() #determine if GAMBLR is running remotely

  if(missing(this_sample_id) & missing(these_samples_metadata)){
    stop("Must provide a single row of metadata or multi-row metadata together with a single sample_id.")
  }else if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata() 
  }
  
  if(nrow(these_samples_metadata) > 1){ 
    #if more than one row was provided, filter to the one we want 
    #just in case more than one row was provided
    if(missing(this_sample_id)){
      stop("Please provide a single sample_id to subset these_samples_metadata to a single row.")
    }else{
      these_samples_metadata = these_samples_metadata %>%
        dplyr::filter(seq_type != "mrna") %>%
        dplyr::filter(sample_id == this_sample_id)
      # If the length is still >1, then >1 seq_type is present for this sample
      if(nrow(these_samples_metadata) > 1){
        stop(glue::glue("There is >1 seq_type for {this_sample_id} in the metadata.
                        Please subset these_samples_metadata to a single seq_type/sample_id combination."))
      }
      if(nrow(these_samples_metadata) == 0){
        stop(glue::glue("No sample with id {this_sample_id} found in the metadata."))
      }
  }}
  
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
    path_template = check_config_and_value(
        "results_flatfiles$ssm$template$clustered$augmented",
        config_name="default")
    path_complete = unname(unlist(glue::glue(path_template)))
    aug_maf_path = paste0(check_config_and_value("project_base",config_name="default"), path_complete)
    local_aug_maf_path = paste0(check_config_and_value("project_base"), path_complete)

  }else{
    warning("Currently the only flavour available to this function is 'clustered'")
  }

  #stop()
  if(remote_session){
    if(verbose){
      if(augmented){
        print(paste("seeking both files:",local_full_maf_path, local_aug_maf_path))
      }else{
        print(paste("seeking file:",local_full_maf_path))
      }
    }
    # Always sync all available files regardless of what the user requested
    # so we don't have to check if they exist remotely every time
    
    if(!file.exists(local_aug_maf_path) & !file.exists(local_full_maf_path)){
      #assume we need to sync everything that's available remotely
      if(verbose){
        print(paste("Missing both files:",local_aug_maf_path,local_full_maf_path))
      }
      #makes an ssh session only when necessary
      remote_session = check_remote_configuration(auto_connect = TRUE)
      #need to obtain a copy of the remote file(s)
      #check if remote file actually exists
      # augmented
      dirN = dirname(local_aug_maf_path)
      status = ssh::ssh_exec_internal(ssh_session,
                                   command=paste("stat",aug_maf_path),
                                   error=F)$status
      if(status==0){
        if(verbose){
              print(paste("found:",aug_maf_path))
              print(paste("local home:",local_aug_maf_path))
        }
        suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
        ssh::scp_download(ssh_session,aug_maf_path,dirN)
      }else{
            #assume augmented MAF doesn't exist remotely and move on
            if(verbose){
              print(paste("not found remotely:",aug_maf_path))
            }
      }
      

      # full
      dirN = dirname(local_full_maf_path)
      status = ssh::ssh_exec_internal(ssh_session,
                                   command=paste("stat",full_maf_path),
                                   error=F)$status
      if(status==0){
          if(verbose){
            print(paste("found:",full_maf_path))
            print(paste("local home:",local_full_maf_path))
        }
        suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
        ssh::scp_download(ssh_session,full_maf_path,dirN)
      }else{
        print(ssh_session)
        stop("failed getting full MAF in get_ssm_by_sample")
      }
    }
    aug_maf_path = local_aug_maf_path
    full_maf_path = local_full_maf_path
    
  }
  if(augmented){
    if(file.exists(aug_maf_path)){
      if(verbose){
        print(paste("setting full_maf_path to",aug_maf_path))
      }
      full_maf_path = aug_maf_path
    }
  }
  # Check if maf exists; if not, return empty data frame
  if(!file.exists(full_maf_path)){
    print(paste("missing: ", full_maf_path))
    message(paste("warning: file does not exist, skipping it.", full_maf_path))
    return()
  }
  sample_ssm = fread_maf(full_maf_path)
  if(min_read_support){
    # drop poorly supported reads but only from augmented MAF
      sample_ssm = dplyr::filter(sample_ssm, t_alt_count >= min_read_support)
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
