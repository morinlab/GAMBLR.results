#' @title Get SSM By Region.
#'
#' @description Retrieve all SSMs from the GAMBL database within a single genomic coordinate range.
#'
#' @details This function lets the user specify a region of interest for returning SSM calls within that region.
#' There are multiple ways a region can be specified. For example, the user can provide the full region in a "region" format (chr:start-end) to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually with `chr`, `start`, and `end` parameters.
#' For more usage examples, refer to the parameter descriptions and examples in the vignettes.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_coding_ssm],
#' [GAMBLR.results::get_coding_ssm_status], [GAMBLR.results::get_ssm_by_patients], [GAMBLR.results::get_ssm_by_sample], [GAMBLR.results::get_ssm_by_samples], [GAMBLR.results::get_ssm_by_regions]
#'
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param qstart Query start coordinate of the range you are restricting to.
#' @param qend Query end coordinate of the range you are restricting to.
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' @param basic_columns Set to FALSE to return MAF with all columns (116). Default is TRUE, which returns the first 45 columns. Note that if streamlined is set to TRUE, only two columns will be returned, regardless of what's specified in this parameter.
#' @param streamlined Return Start_Position and Tumor_Smaple_Barcode as the only two MAF columns. Default is FALSE. Setting to TRUE will overwrite anything specified with `basic_columns`.
#' @param maf_data An already loaded MAF like object to subset to regions of interest.
#' @param this_seq_type The seq_type you want back, default is genome.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#' @param verbose Boolean parameter set to FALSE per default.
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#'
#' @rawNamespace import(vroom, except = c("col_skip", "fwf_positions", "default_locale", "date_names_lang", "cols_only", "output_column", "col_character", "col_guess", "spec", "as.col_spec", "fwf_cols", "cols", "col_date", "col_datetime", "locale", "col_time", "cols_condense", "col_logical", "col_number", "col_integer", "col_factor", "fwf_widths", "date_names_langs", "problems", "date_names", "col_double", "fwf_empty"))
#' @import dplyr RMariaDB DBI stringr glue GAMBLR.helpers
#'
#' @examples
#' #basic usage
#' my_mutations = GAMBLR.results:::get_ssm_by_region(region = "chr8:128723128-128774067")
#'
#'
#' #keep all 116 columns in the read MAF
#' bcl2_all_details = get_ssm_by_region(region = "chr18:60796500-60988073",
#'                                      basic_columns = FALSE)
#'
get_ssm_by_region = function(chromosome,
                             qstart,
                             qend,
                             region = "",
                             these_sample_ids = NULL,
                             these_samples_metadata = NULL,
                             basic_columns = TRUE,
                             streamlined = FALSE,
                             maf_data,
                             this_seq_type = "genome",
                             projection = "grch37",
                             from_indexed_flatfile = TRUE,
                             augmented = TRUE,
                             min_read_support = 3,
                             mode = "slms-3",
                             verbose = FALSE){
  
  
  seq_type = this_seq_type
  if(mode == "strelka2"){
    message("Mode is set to strelka2. Streamlined = TRUE is hardcoded for this mode...")
    streamlined = TRUE #force streamlined to TRUE, if strelka2 output is requested.
    augmented = FALSE #force augmented to FALSE (since t_alt_count column is not available for the strelka2 bed file).
    maf_columns = c("Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode")
    maf_column_types = "iiic"

    #add some checks
    if(projection == "hg38"){
      stop("Strelka2 outputs are currently only available in respect to grch37...")
    }
    if(this_seq_type == "capture"){
      stop("Genome is currently the only available seq_type for strelka2 outputs...")
    }
  }else if(mode == "slms-3"){
    if(streamlined){
      maf_columns = names(maf_header)[c(6, 16, 42)]
      maf_column_types = "ici"

    }else if(basic_columns){ #get first 45 columns of the MAF
      maf_columns = names(maf_header)[c(1:45)]
      maf_column_types =  "ciccciiccccccclcccclllllllllllllllccccciiiiii"
    }else{
      maf_columns = names(maf_header) #return all MAF columns (116)
      maf_column_types = "ciccciiccccccclcccclllllllllllllllccccciiiiiiccccccccccccinnccccccccccccccccccclcccccccccnclcncccclncccclllllllllicn"
    }
  }

  #check that maf_columns requested all exist in the header and get their indexes
  if(!all(maf_columns %in% names(maf_header))){
    stop("Cannot find one of the requested maf_columns in your MAF header")
  }

  #get MAF column indexes
  maf_indexes = maf_header[maf_columns]
  maf_indexes = maf_indexes[order(maf_indexes)]
  maf_columns = names(maf_indexes)
  maf_indexes = unname(maf_indexes)

  #get config values
  tabix_bin = check_config_and_value("dependencies$tabix")
  table_name = check_config_and_value("results_tables$ssm")
  db = check_config_and_value("database_name")
  base_path = check_config_and_value("project_base")
  base_path_remote = check_config_and_value("project_base",config_name="default")

  #get absolute file paths based on the selected mode and check existance for the file
  if(from_indexed_flatfile){
    if(mode == "slms-3"){
      if(augmented){
        maf_partial_path = check_config_and_value("results_flatfiles$ssm$template$merged$augmented")
      }else{
        maf_partial_path = check_config_and_value("results_flatfiles$ssm$template$merged$deblacklisted")
      }
    }else if (mode == "strelka2"){
      maf_partial_path = check_config_and_value("results_flatfiles$ssm$all$strelka2")
    }else{
      stop("You requested results from indexed flatfile. The mode should be set to either slms-3 (default) or strelka2. Please specify one of these modes.")
    }

    #use glue to get the absolute path
    maf_path = glue::glue(maf_partial_path)
    full_maf_path = paste0(base_path, maf_path)
    
    if(mode == "slms-3"){
      full_maf_path_comp = paste0(full_maf_path, ".bgz")

    }else if(mode == "strelka2"){
      full_maf_path_comp = gsub('.{3}$', 'bed', full_maf_path) #do we instead want to add the exact path to the file in the config, or is this acceptable?
      full_maf_path_comp = paste0(full_maf_path_comp, ".gz")
    }

    #check if file is existing or missing
    if(!file.exists(full_maf_path_comp)){
      print(paste("missing:", full_maf_path_comp))
      check_host(verbose=TRUE)
      if(verbose){
        print("using local file")
        print(paste("HERE:",full_maf_path_comp))
      }
      stop("failed to find the file needed for this")
    }
  }

  #split region into chunks (chr, start, end) and deal with chr prefixes based on the selected projection
  if(!region == ""){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    if(projection == "grch37"){
      region = stringr::str_replace(region, "chr", "")
    }
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = as.numeric(startend[1])
    qend = as.numeric(startend[2])
  }else{
    if(projection =="grch37"){
      chromosome = gsub("chr", "", chromosome)
    }
    region = paste0(chromosome, ":", qstart, "-", qend)
  }

  if(projection =="grch37"){
    chromosome = gsub("chr", "", chromosome)
  }

  #use vroom on indexed maf files (if maf_data is not provided)
  if(missing(maf_data)){
    
    #get tabix command
    if(mode == "slms-3"){
      tabix_command = paste(tabix_bin, full_maf_path_comp, region, "| cut -f", paste(maf_indexes, collapse = ","))
    }else if(mode == "strelka2"){
      tabix_command = paste(tabix_bin, full_maf_path_comp, region)
    }
    if(verbose){
        print(tabix_command)
    }

    #execute the tabix command
    muts = system(tabix_command, intern = TRUE)
    if(verbose){
      print(paste("TYPES:"))
      print(maf_column_types)
      print("NAMES:")
      print(maf_columns)
      print("NUM:")
      print(length(muts))
    }

    if(length(muts)==0){
      maf_types_sep = str_split(maf_column_types, pattern = "")[[1]] %>%
          str_replace_all("c", "character") %>%
          str_replace_all("l|i|n", "numeric")
      if(verbose){
        print("adding dummy line")
      }

      muts_region = read.table(textConnection(""), col.names = maf_columns, colClasses = maf_types_sep)

      
    }else{
      #this is what gets executed with default parameters
      muts_region = vroom::vroom(I(muts), col_types = paste(maf_column_types, collapse = ""), col_names = maf_columns, delim = "\t")
    }

    if(verbose){
      print('SUCCESS')
    }
    if(augmented){
        # drop poorly supported reads but only from augmented MAF
      muts_region = dplyr::filter(muts_region, t_alt_count >= min_read_support)
    }
  
  }else{
    muts_region = dplyr::filter(maf_data, Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
    muts_region = dplyr::filter(maf_data, Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
  }
  if(!missing(these_sample_ids)){
    stop("deprecated argument: these_sample_ids, please use these_samples_metadata instead")
  }
  if(!missing(these_samples_metadata)){
    sample_ids = pull(these_samples_metadata,sample_id)
    # subset sample ids
    muts_region = dplyr::filter(muts_region, Tumor_Sample_Barcode %in% sample_ids)
  }
  if(streamlined){
    muts_region = muts_region %>%
      dplyr::select(Start_Position, Tumor_Sample_Barcode)
  }

  return(muts_region)
}
