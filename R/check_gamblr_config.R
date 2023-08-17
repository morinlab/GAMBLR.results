#' @title Check GAMBLR Config.
#'
#' @description Check that the GAMBLR config you have loaded will work in your setup.
#'
#' @details This function is mostly for remote GAMBLRs to ensure they keep their local mirror of the GAMBL data up-to-date.
#'
#' @param compare_timestamps Whether the function will compare timestamps on your local files to the remote copy. Only relevant if you are working remotely.
#' @param ssh_session The ssh_session object see get_ssh_session() for more information. Only relevant if you are working remotely.
#' @param archive_mode This is not currently working but the idea here is to keep a GSC archive of GAMBL in sync with the actively updated outputs
#' @param force_backup Boolean parameter set to FALSE per default.
#'
#' @import dplyr tidyr stringi glue
#' @export
#'
#' @examples
#' check_gamblr_config()
#'
check_gamblr_config = function(compare_timestamps=FALSE,
                               ssh_session,
                               archive_mode=FALSE,
                               force_backup=FALSE){
  files_to_check = c()
  #get all the wildcards we'll need
  seq_type = get_template_wildcards("seq_types")
  seq_type_filter = seq_type
  unix_group = get_template_wildcards("unix_groups")
  projection = get_template_wildcards("projections")

  #data frame for seq_type/projection expansion
  projection_expanded = tidyr::expand_grid(seq_type = get_template_wildcards("seq_types"),projection = get_template_wildcards("projections"))
  print(projection_expanded)
  #resources section of config (only needs blacklist right now)
  blacklist_f = check_config_value(config::get("resources")$blacklist$template)
  blacklist_f = mutate(projection_expanded,output=glue::glue(blacklist_f)) %>% pull(output)
  files_to_check = c(files_to_check,blacklist_f)

  merged_keys = names(check_config_value(config::get("results_merged")))
  #skip any file starting with "/"
  for (merge in merged_keys){
    merge_path = config::get("results_merged")[merge]
    if(!grepl("^/",merge_path)){
      files = unlist(merge_path)
      #print(names(files))

      for(f in files){
        #print(paste("CHECKING",f))
        if(grepl("\\{",f)){
          #print("contains wildcards, using all wildcards:")
          if(stringi::stri_count_fixed(f,"{")>1){
            print("Multiple wildcards!")
            wildcards = grob_wildcards(f)
            print(wildcards)
            if("projection" %in% wildcards & "seq_type" %in% wildcards){
              #use the same expansion approach we used above
              print(f)
              f = mutate(projection_expanded,output=glue::glue(f)) %>% pull(output)
              files_to_check = c(files_to_check,f)
            }else{
              warning("Unrecognized wildcards. Not sure how to handle this, using default approach")
              flavour =get_template_wildcards("results_merged",names(files))
              all_f = glue::glue(f)
              print(all_f)
              files_to_check = c(files_to_check,all_f)
            }
          }else{
            flavour =get_template_wildcards("results_merged",names(files))
            all_f = glue::glue(f)
            #print(all_f)
            files_to_check = c(files_to_check,all_f)
          }
        }else{
          files_to_check = c(files_to_check,f)
        }

      }

    }
  }
  print(f)
  mia=check_file_details(files_to_check)
  l_missing = length(mia)
  if(l_missing){
    warning(paste("There were",l_missing,"files that cannot be found (see above). If this is unexpected, try to obtain them."))
    print("MISSING FILES:")
    print(mia)
  }
  if(compare_timestamps){
    check_times(files_to_check,archive_mode,force_backup)
  }
  print("DONE!")
}
