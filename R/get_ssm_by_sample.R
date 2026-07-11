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
#' only slms_3 is supported).
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
#' @param variant_classification_filter Optional character vector of
#' Variant_Classification values to keep (e.g. GAMBLR.helpers::vc_nonSynonymous).
#' When supplied, the MAF is pre-filtered with grep before it is ever parsed
#' by R: locally, this means non-matching lines are never type-converted
#' into a data frame; when the file must be fetched from a remote host over
#' the existing ssh_session, the grep is instead run on the remote host so
#' only the (pre-filtered) matching lines are transferred over the network
#' rather than the whole MAF. This avoids the memory, CPU and (for remote
#' sessions) network cost of materializing the full MAF when, as is typical,
#' only a small fraction of variants match. grep matches whole words
#' anywhere on the line, so it is a fast superset filter; an exact
#' dplyr::filter on Variant_Classification is always applied afterwards to
#' guarantee correctness. Default NULL returns all variant classifications
#' (original behaviour, and no ssh optimization is attempted).
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
#'  # only load non-synonymous variants, avoiding the cost of materializing
#'  # the full MAF when most rows will just be discarded
#'  nonsyn_maf = GAMBLR.results:::get_ssm_by_sample(
#'    get_gambl_metadata() %>% dplyr::filter(sample_id=="13-27975_tumorA"),
#'    variant_classification_filter = GAMBLR.helpers::vc_nonSynonymous
#'  )
#'
get_ssm_by_sample = function(these_samples_metadata,
                             tool_name = "slms_3",
                             projection = "grch37",
                             augmented = TRUE,
                             flavour = "clustered",
                             min_read_support = 3,
                             basic_columns = TRUE,
                             maf_cols = NULL,
                             variant_classification_filter = NULL,
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

  }else if(flavour == "sage"){
    vcf_base_name = "sage.combined"
    path_template = check_config_and_value(
      "results_flatfiles$ssm$template$sage$deblacklisted",
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
        "results_flatfiles$ssm$template$sage$augmented",
        config_name="default")
    path_complete = unname(unlist(glue::glue(path_template)))
    aug_maf_path = paste0(check_config_and_value("project_base",config_name="default"), path_complete)
    local_aug_maf_path = paste0(check_config_and_value("project_base"), path_complete)

  }
  else{
    warning("Currently the only flavour available to this function is 'clustered'")
  }

  #stop()
  # When variant_classification_filter is set, build a grep command that
  # cheaply drops lines that can't match (grep matches whole words anywhere
  # on the line, so this is a fast superset; the exact dplyr::filter below
  # is what guarantees correctness either way). type="sh" keeps the quoting
  # POSIX-style regardless of the local platform, since this command may be
  # run on a remote host over ssh.
  grep_classification_cmd = function(path, classifications){
    terms = paste(shQuote(classifications, type = "sh"), collapse = " -e ")
    paste("grep -h -w -F -e Hugo_Symbol -e", terms, shQuote(path, type = "sh"))
  }
  sample_ssm_loaded = FALSE
  
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

      if(!is.null(variant_classification_filter)){
        # Only the single MAF that will actually be used is needed, and
        # only the pre-filtered matching lines need to cross the network -
        # run grep on the remote host itself instead of scp-ing the whole
        # (unfiltered) file down first.
        remote_path = full_maf_path
        if(augmented){
          aug_status = ssh::ssh_exec_internal(ssh_session,
                                       command=paste("stat",aug_maf_path),
                                       error=F)$status
          if(aug_status==0) remote_path = aug_maf_path
        }
        print(paste("grep:",remote_path))
        cmd = grep_classification_cmd(remote_path, variant_classification_filter)
        
        if(verbose){
          print(paste("running on remote host:", cmd))
        }
        result = ssh::ssh_exec_internal(ssh_session, command = cmd, error = FALSE)
        # grep exits 1 (not an error) when nothing matches; anything else
        # (e.g. file not found on the remote host) is a real failure.
        if(!result$status %in% c(0,1)){
          stop(paste("failed getting filtered MAF in get_ssm_by_sample from", remote_path))
        }
        sample_ssm = fread(text = rawToChar(result$stdout))
        sample_ssm_loaded = TRUE
      }else{
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
  if(!sample_ssm_loaded){
    # Check if maf exists; if not, return empty data frame
    if(!file.exists(full_maf_path)){
      message(paste("warning: file does not exist, skipping it.", full_maf_path))
      return()
    }
    if(!is.null(variant_classification_filter)){
      sample_ssm = fread(cmd = grep_classification_cmd(full_maf_path, variant_classification_filter))
    }else{
      # fread_maf() applies a rigid, hand-maintained colClasses spec via
      # readr that doesn't necessarily match every flavour's actual column
      # contents (e.g. it mis-types some sage MAF columns), which can
      # silently drop every row. We already have the header in the file
      # itself, so just let fread infer types directly - consistent with
      # the grep-filtered path above.
      sample_ssm = fread(full_maf_path)
    }
  }
  # fread() infers column types per-file; for sparsely-populated columns
  # (e.g. gnomAD_* frequencies, usually NA) this can disagree across samples
  # (logical in one file, character/double in another), which breaks
  # bind_rows() when per-sample MAFs are later combined (e.g. by
  # get_ssm_by_samples()). Coerce every column present to its canonical type
  # (GAMBLR.utils::maf_column_classes(), the same map fread_maf() uses) after
  # the lenient read, rather than constraining the read itself -- this can
  # only turn a genuinely non-conforming value into NA, it never drops rows
  # the way applying colClasses at read time did.
  col_types = GAMBLR.utils::maf_column_classes()
  for(col in intersect(names(col_types), names(sample_ssm))){
    target = col_types[[col]]
    sample_ssm[[col]] = suppressWarnings(switch(target,
      character = as.character(sample_ssm[[col]]),
      integer   = as.integer(sample_ssm[[col]]),
      numeric   = as.numeric(sample_ssm[[col]]),
      logical   = as.logical(sample_ssm[[col]]),
      sample_ssm[[col]]
    ))
  }
  if(!is.null(variant_classification_filter)){
    # grep above is a superset match (whole word, anywhere on the line);
    # this exact filter is what guarantees correctness.
    sample_ssm = dplyr::filter(sample_ssm, Variant_Classification %in% variant_classification_filter)
  }
  if(min_read_support){
    # drop poorly supported reads, but never drop a row just because
    # t_alt_count itself is missing (e.g. some flavours/deblacklisted
    # MAFs don't populate read-depth columns at all)
      sample_ssm = dplyr::filter(sample_ssm, t_alt_count >= min_read_support | is.na(t_alt_count))
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
