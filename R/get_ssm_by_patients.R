#' @title Get SSM By Patients.
#'
#' @description Get MAF-format data frame for more than one patient.
#'
#' @details This function returns variants from a set of patients avoiding duplicated mutations from multiple samples from that patient (i.e. unique superset of variants).
#' This is done either by combining the contents of individual MAF files or subset from a merged MAF (wraps [GAMBLR.results::get_ssm_by_samples]).
#' In most situations, this should never need to be run with `subset_from_merge = TRUE`. Instead use one of [GAMBLR.results::get_coding_ssm] or [GAMBLR.results::get_ssm_by_region].
#' This function expects either a vector of patient IDs (`thse_patients_ids`) or an already subset metadata table (`these_samples_metadata`).
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_coding_ssm], [GAMBLR.results::get_coding_ssm_status],
#' [GAMBLR.results::get_ssm_by_sample], [GAMBLR.results::get_ssm_by_samples], [GAMBLR.results::get_ssm_by_region], [GAMBLR.results::get_ssm_by_regions]
#'
#' @param these_patient_ids A vector of patient IDs that you want results for. The user can also use a metadata table that has been subset to the patient IDs of interest (`these_samples_metadata`).
#' @param these_samples_metadata A metadata subset to contain the rows corresponding to the patients of interest. If the vector of patient IDs is missing (`these_patient_ids`), this function will default to all patient IDs in the metadata table given to this parameter.
#' @param tool_name Only supports slms-3 currently.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param seq_type The seq type you want results for. Default is "genome".
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatibility with additional variant calling parameters/versions.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters (matching columns in MAF).
#' @param subset_from_merge Instead of merging individual MAFs, the data will be subset from a pre-merged MAF of samples with the specified `seq_type`.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead.
#' @param engine Specify one of readr or fread_maf (default) to change how the large files are loaded prior to subsetting. You may have better performance with one or the other but for me fread_maf is faster and uses a lot less RAM.
#'
#' @return A data frame with SSM calls for the selected patients in MAF format.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' #example 1, using a vector of patient IDs.
#' patients = c("00-14595", "00-15201", "01-12047")
#' patients_maf = get_ssm_by_patients(these_patient_ids = patients,
#'                                    seq_type = "genome",
#'                                    subset_from_merge = FALSE)
#'
#' #example 2, using a metadata table, subset to the patient IDs of interest.
#' patient_meta = get_gambl_metadata(seq_type_filter = "genome")
#' patient_meta = dplyr::filter(patient_meta, patient_id %in% patients)
#'
#' patients_maf_2 = get_ssm_by_patients(these_samples_metadata = patient_meta,
#'                                      subset_from_merge = FALSE)
#'
get_ssm_by_patients = function(these_patient_ids,
                               these_samples_metadata,
                               tool_name = "slms-3",
                               projection = "grch37",
                               seq_type = "genome",
                               flavour = "clustered",
                               min_read_support = 3,
                               basic_columns = TRUE,
                               maf_cols = NULL,
                               subset_from_merge = FALSE,
                               augmented = TRUE,
                               engine='fread_maf'){

  check_remote_configuration(auto_connect = TRUE)
  if(!subset_from_merge){
    message("WARNING: on-the-fly merges can be extremely slow and consume a lot of memory. Use at your own risk. ")
  }
  if(missing(these_patient_ids)){
    if(missing(these_samples_metadata)){
      stop("must supply either a vector of patient_ids or the metadata for those patients as these_samples_metadata")
    }
    these_patient_ids = pull(these_samples_metadata,patient_id) %>% unique()
  }
  augmented = TRUE
  #always requires augmented MAFs to ensure all variants from the patient are included
  to_exclude = get_excluded_samples(tool_name)
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = seq_type) %>%
      dplyr::filter(patient_id %in% these_patient_ids) %>%
      dplyr::filter(!sample_id %in% to_exclude)
  }else{
    these_samples_metadata = these_samples_metadata %>%
      dplyr::filter(patient_id %in% these_patient_ids) %>%
      dplyr::filter(!sample_id %in% to_exclude)
  }
  # Removed because this drops all but one sample for the patient!
  #these_samples_metadata = group_by(these_samples_metadata, patient_id) %>%
  #  slice_head() %>%
  #  ungroup()

  these_sample_ids = pull(these_samples_metadata, sample_id)

  return(get_ssm_by_samples(these_sample_ids = these_sample_ids,
                            these_samples_metadata = these_samples_metadata,
                            tool_name = tool_name,
                            projection = projection,
                            seq_type = seq_type,
                            flavour = flavour,
                            min_read_support = min_read_support,
                            subset_from_merge = subset_from_merge,
                            augmented = augmented,
                            basic_columns = basic_columns,
                            maf_cols = maf_cols,
                            engine=engine))
}
