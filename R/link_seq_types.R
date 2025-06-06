#' @title Link seq type IDs
#' 
#' @description Link sample IDs of a seq type to sample IDs of a different seq type that share a same
#'   feature, *e.g.* same `biopsy_id`.
#' 
#' @details If there is more than one sample id of the desired seq type associated with a same value
#'  in `intermediary`, the `multi_desired_id` column of the output table will show 1 instead of 0.
#' 
#' @param these_sample_ids A vector of sample IDs of the same seq type specified in `given_seq_type`.
#' @param given_seq_type A single string specifying the seq type of `these_sample_ids`. Possible 
#'   values are one of "genome", "capture", or "mrna".
#' @param desired_seq_type A single string specifying the seq type of the desired sample IDs.
#'   Possible values are one of "genome", "capture", "mrna", or "any". If "any", sample IDs
#'   of any other seq type that is not `given_seq_type` will be returned. 
#' @param intermediary A vector of one or more strings with column names in the metadata which will 
#'   act as intermediaries to link the given sample IDs with the desired sample IDs.
#'
#' @return A data frame with the given and desired sample IDs, values of the feature specified in
#'   the `intermediary` parameter, and the `multi_desired_id` column.
#' 
#' @import dplyr
#' @export
#'
#' @examples 
#' \dontrun{
#' link_seq_types(
#'   these_sample_ids = c("00-14595_tumorA", "04-21622_tumorB", "01-20774T"), 
#'   given_seq_type = "genome",
#'   desired_seq_type = "mrna", 
#'   intermediary = c("patient_id", "biopsy_id")
#' )
#' }
#' @keywords internal
link_seq_types <- function(these_sample_ids, 
                           given_seq_type, 
                           desired_seq_type, 
                           intermediary){
  
  # check parameters
  stopifnot("`given_seq_type` must be one of \"genome\", \"capture\", or \"mrna\"." = 
              length(given_seq_type) == 1 & given_seq_type %in% c("genome", "capture", "mrna"))
  stopifnot("`desired_seq_type` must be one of \"genome\", \"capture\", \"mrna\" or \"any\"." = 
              length(desired_seq_type) == 1 & desired_seq_type %in% c("genome", "capture", "mrna"))
  stopifnot("given_seq_type and desired_seq_type must be different." = 
              given_seq_type != desired_seq_type)
  
  # get metadata
  meta = get_gambl_metadata(seq_type_filter = c(given_seq_type, desired_seq_type))
  missing_columns <- ! intermediary %in% names(meta)
  if(any(missing_columns)){
    k <- intermediary[missing_columns] %>% 
      paste(collapse = ", ") %>% 
      gettextf("`intermediary` must be column names present in the metadata. Missing coumns: %s.", .)
    stop(k)
  }
  
  # subset metadata
  meta = dplyr::select(meta, sample_id, seq_type, all_of(intermediary)) %>% 
    dplyr::filter(seq_type %in% c(given_seq_type, desired_seq_type))
  
  # add intermediary columns to the input these_sample_ids
  these_sample_ids = dplyr::filter(meta, seq_type == given_seq_type) %>% 
    dplyr::select(sample_id, all_of(intermediary)) %>% 
    left_join( as.data.frame(these_sample_ids), .,
               by = join_by(these_sample_ids == sample_id) )
  
  # add the desired ids to these_sample_ids
  these_sample_ids = dplyr::filter(meta, seq_type == desired_seq_type) %>% 
    dplyr::select(-seq_type) %>% 
    unique %>% 
    left_join(these_sample_ids, ., by = intermediary, relationship = "many-to-many")
  
  # rename columns
  given_seq_type_column_name <- paste0(given_seq_type, "_sample_id")
  desired_seq_type_column_name <- paste0(desired_seq_type, "_sample_id")
  these_sample_ids = rename( these_sample_ids, {{given_seq_type_column_name}} := these_sample_ids,
          {{desired_seq_type_column_name}} := sample_id )
  
  # add multi_desired_id column
  these_sample_ids = split(these_sample_ids[,desired_seq_type_column_name], 
                           these_sample_ids[,given_seq_type_column_name]) %>%
    .[lengths(.) > 1] %>%
    names %>%
    { mutate(these_sample_ids, 
             multi_desired_id = ifelse(these_sample_ids[[given_seq_type_column_name]] %in% ., 1, 0)) }
  
  these_sample_ids
}
