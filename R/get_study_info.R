#' @title Get Study Info.
#'
#' @description Function for retrieving study specific identifiers.
#'
#' @details This function takes one required parameter (`dir`).
#' This is the relative path to the main directory for the study of interest.
#' The function reads in the meta_study.txt file and extracts unique study identifiers.
#' This returns a list that holds all the identifiers.
#' The user can also return the study identifiers to the global environment.
#' To do so, set the `list_to_global = TRUE`.
#'
#' @param dir The relative path to the study directory (expects to find meta_study.txt in this folder).
#' @param list_to_global Boolean parameter, if set to TRUE all study identifiers will be returned to the global environment. Default is FALSE.
#'
#' @return A list with study specific identifiers, or nothing (i.e list_to_global = TRUE).
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' #return study identifiers as  list:
#' my_study_info= get_study_info(dir = "path/to/study/")
#'
#' #return all identifiers to the global environment:
#' get_study_info(dir = "path/to/study/", list_to_global = TRUE)
#' }
#' 
get_study_info = function(dir,
                          list_to_global = FALSE){

  #check if meta data file exists in the selected directory
  if(file.exists(paste0(dir, "meta_study.txt"))){
    study_meta =  data.table::fread(file = paste0(dir, "meta_study.txt"), sep = "\t", header = FALSE)
  }else{
    stop("Unable to find meta_study.txt in the specified folder (dir)...")
  }

  #tidy the data frame
  meta_info = gsub(".*: ","", study_meta$V1)

  #create a list with study identifiers
  meta_list = list(cancer_type = meta_info[1],
                   project_name = meta_info[2],
                   human_friendly_name = meta_info[3],
                   short_name = meta_info[4],
                   description = meta_info[5])

  #add vectors to global environment
  if(list_to_global){
    list2env(meta_list, envir = .GlobalEnv)
  }else{
    return(meta_list)
  }
}
