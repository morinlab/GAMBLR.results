#' @title Load MutationTimeR Results
#'
#' @description This function retrieves and loads the MutationTimeR
#' timing estimates for a single sample_id for all successfully timed
#' and un-timed SSMs and CNAs.
#'
#' @param this_sample_metadata Metadata with one row
#' containing the details for the desired genome sample_id
#' @param projection Genome build projection (e.g., "hg38" or "grch37").
#' @param verbose Set to TRUE for a chattier experience. Default is FALSE.
#' @import readr dplyr
#' @export
#' @return a named list containing two data.frames with
#' the ssm and cna timing information
get_timed_mutations <- function(this_sample_metadata, projection, verbose = FALSE) {
  if (nrow(this_sample_metadata) > 1) {
    stop("this_sample_metadata must contain exactly one row")
  }
  if (!projection %in% c("hg38", "grch37")) {
    stop("Genome build projection must be either hg38 or grch37")
  }
  remote_session = NULL
  normal_sample_id <- pull(this_sample_metadata, normal_sample_id)
  sample_id <- pull(this_sample_metadata, sample_id)
  unix_group <- pull(this_sample_metadata, unix_group)
  tumour_sample_id <- sample_id
  pattern_cna <- "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/{unix_group}/mutationtimer-1.0/99-outputs/timed_cna/{projection}/{sample_id}--{normal_sample_id}_timed_cna.{projection}.tsv"
  input_cna <- glue::glue(pattern_cna)
  pattern_ssm <- "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/{unix_group}/mutationtimer-1.0/99-outputs/timed_ssm/{projection}/{sample_id}--{normal_sample_id}_timed_ssm.{projection}.tsv"
  input_ssm <- glue::glue(pattern_ssm)
  local_base = config::get()$project_base
  local_cna <- "{local_base}/{unix_group}/mutationtimer-1.0/99-outputs/timed_cna/{projection}/{sample_id}--{normal_sample_id}_timed_cna.{projection}.tsv"
  local_ssm <- "{local_base}/{unix_group}/mutationtimer-1.0/99-outputs/timed_ssm/{projection}/{sample_id}--{normal_sample_id}_timed_ssm.{projection}.tsv"
  local_cna = glue::glue(local_cna)
  local_ssm = glue::glue(local_ssm)
  if(!file.exists(input_cna) | !file.exists(input_ssm)){
    #might be a remote session
    if(file.exists(local_cna) & file.exists(local_ssm)){
        if(verbose) {
          print("found both files locally!")
        }
        input_cna = local_cna
        input_ssm = local_ssm
    }else{
        #try to get the files
        remote_session <- check_remote_configuration(auto_connect = TRUE)
    }
  }

  if (!is.null(remote_session)) {
    # check if file exists
    status_cna <- ssh::ssh_exec_internal(ssh_session,
                                            command = paste("stat", input_cna),
                                            error = F)$status
    status_ssm <- ssh::ssh_exec_internal(ssh_session,
                                             command = paste("stat", input_ssm),
                                            error = F)$status
    # first check if we already have a local copy
    # Load data from local copy or get a local copy from the remote path first
    if (status_cna == 0 & status_ssm == 0) {
      if (verbose) {
        print(paste("found:", input_cna))
        print(paste("local home:", local_cna))
      }
      dirN <- dirname(local_cna)

      suppressMessages(suppressWarnings(dir.create(dirN, recursive = T)))
      if (!file.exists(local_cna)) {
        if(verbose) {
          print(paste("DOWNLOADING:", input_cna))
        }
        ssh::scp_download(ssh_session, input_cna, dirN)
      }
      dirN <- dirname(local_ssm)
      suppressMessages(suppressWarnings(dir.create(dirN, recursive = T)))
      if (!file.exists(input_ssm)) {
        if(verbose) {
          print(paste("DOWNLOADING:", input_ssm))
        }
        ssh::scp_download(ssh_session, input_ssm, dirN)
      }
      input_ssm = local_ssm
      input_cna = local_cna
    }
  }
  cna = suppressMessages(read_tsv(input_cna, progress = FALSE))
  cna = create_genomic_data(cna,projection)
  ssm = suppressMessages(read_tsv(input_ssm, progress = FALSE))
  ssm = create_genomic_data(ssm,projection)
  return(list(CNA=cna,SSM=ssm))
}
