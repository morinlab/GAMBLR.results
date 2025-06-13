#' @title Get SSM By Region.
#'
#' @description Retrieve all SSMs from the GAMBL database within a single genomic coordinate range.
#'
#' @details This function lets the user specify a region of interest
#' for obtaining SSM calls within that region. In most cases, you
#' should be using [GAMBLR.results::get_ssm_by_regions] instead.
#' There are multiple ways a region can be specified. For example,
#' the user can provide the full region in a "region" format (chr:start-end)
#' to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually
#' with `chr`, `start`, and `end` parameters.
#'
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param qstart Query start coordinate of the region of interest.
#' @param qend Query end coordinate of the region of interest.
#' @param region Full region definition specified as a character vector e.g. "chrX:1234-5678"
#'  instead of specifying chromosome, start and end separately.
#' @param these_samples_metadata Optional metadata table.
#'  If provided, it will return SSM calls for the only the samples in the metadata table.
#'  Otherwise it will use all samples from `get_gambl_metadata()` of the appropriate seq_types.
#' @param maf_data Use an already loaded MAF data frame.
#'  If you would like all columns of this input maf returned,
#'  set `streamlined = FALSE` and `basic_columns = FALSE`.
#'  Otherwise the first 45 columns will be returned.
#' @param streamlined If TRUE, only 3 columns will be returned:
#'  start, sample_id, and region in the format "chr:start-end". Default is FALSE.
#'  Note: if this parameter is TRUE, the function will disregard anything specified with `basic_columns`.
#' @param basic_columns Parameter to be used when streamlined is FALSE.
#'  Set this parameter to TRUE (default) to return a MAF with the standard 45 columns.
#'  Set to FALSE to return a MAF with all columns (116).
#'  If you provided `maf_data` with more than 45 columns, set to FALSE to return all columns of
#'  `maf_data`, otherwise it will return the first 45.
#' @param tool_name Accepts either "slms_3" (default) or "strelka2"
#'  (forces `streamlined=TRUE`) to indicate which variant caller to use. Note: strelka2 will force
#'  `augmented=FALSE` as data is not available in that case. Formerly called "mode".
#' @param augmented Default: TRUE. Setting to FALSE will subset from a pre-merged MAF.
#'  Obtaining variants from the original MAFs for each sample would not be computationally efficient here.
#'  Provide to `maf_data` the results of `get_ssm_by_samples` if those are needed.
#' @param projection Obtain variants projected to this reference, one of grch37 (default) or hg38.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count
#'  (for cleaning up augmented MAFs). Default: 3.
#' @param verbose Boolean parameter set to FALSE per default.
#' @param this_seq_type Deprecated. Inferred from these_samples_metadata
#' @param these_sample_ids Deprecated. Inferred from these_samples_metadata
#'
#'
#' @return A data frame of variants in 3 column format or in MAF-like format (one row per mutation).
#'
#' @rawNamespace import(vroom, except = c("col_skip", "fwf_positions", "default_locale", "date_names_lang", "cols_only", "output_column", "col_character", "col_guess", "spec", "as.col_spec", "fwf_cols", "cols", "col_date", "col_datetime", "locale", "col_time", "cols_condense", "col_logical", "col_number", "col_integer", "col_factor", "fwf_widths", "date_names_langs", "problems", "date_names", "col_double", "fwf_empty"))
#' @import dplyr stringr glue GAMBLR.helpers
#'
#' @examples
#' #basic usage
#' genomes = get_gambl_metadata() %>%
#'              dplyr::filter(seq_type %in% "genome")
#' my_mutations = GAMBLR.results:::get_ssm_by_region(
#'             region = "chr8:128723128-128774067",
#'             these_samples_metadata = genomes)
#'
#'
#' #keep all 116 columns in the read MAF
#' bcl2_all_details = GAMBLR.results:::get_ssm_by_region(
#'             region = "chr18:60796500-60988073",
#'             these_samples_metadata = genomes,
#'             basic_columns = FALSE)
#' @keywords internal
get_ssm_by_region = function(chromosome,
                             qstart,
                             qend,
                             region = "",
                             these_samples_metadata,
                             maf_data,
                             streamlined = FALSE,
                             basic_columns = TRUE,
                             tool_name = "slms_3",
                             augmented = TRUE,
                             projection = "grch37",
                             min_read_support = 3,
                             verbose = FALSE,
                             these_sample_ids,
                             this_seq_type){

  if(!missing(this_seq_type) | !missing(these_sample_ids)){
    stop("this_seq_type and these_sample_ids are deprecated. Use these_samples_metadata instead")
  }

  if(!projection %in% c("grch37", "hg38")){
    stop("projection must be either grch37 or hg38")
  }

  if(length(tool_name) != 1){
    stop("tool_name can only be a single value, either slms_3 or strelka2")
  }else if(!tool_name %in% c("slms_3", "strelka2")){
    stop("tool_name must be either slms_3 or strelka2")
  }

  to_exclude = get_excluded_samples(tool_name)

  if(missing(these_samples_metadata)){
    # kept for legacy, assumes user provided these_sample_ids
    message("CAUTION! these_samples_metadata was not provided. Using all of get_gambl_metadata().")
    these_samples_metadata = get_gambl_metadata() %>%
      dplyr::filter(seq_type!="mrna") %>%
      dplyr::filter(!sample_id %in% to_exclude)

  }else{
    #drop unsupported seq_type and samples to exclude
    these_samples_metadata = dplyr::filter(these_samples_metadata,seq_type!="mrna") %>%
        dplyr::filter(!sample_id %in% to_exclude)
  }

  if(tool_name == "strelka2"){
    message("tool_name is set to strelka2. Streamlined = TRUE and augmented = FALSE are hardcoded for this tool...")
    streamlined = TRUE # force streamlined to TRUE, if strelka2 output is requested.
    augmented = FALSE # force augmented to FALSE (since t_alt_count column is not available for the strelka2 bed file).
    maf_columns = c("Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode")
    maf_column_types = "iiic"

    #add some checks
    if(projection == "hg38"){
      stop("Strelka2 outputs are currently only available in respect to grch37...")
    }
    if(length(unique(these_samples_metadata$seq_type))>1){
      stop("More than one seq_type in these_samples_metadata. Genome is currently the only available seq_type for strelka2 outputs...")
    }else if(unique(these_samples_metadata$seq_type) != "genome"){
      stop("Genome is currently the only available seq_type for strelka2 outputs...")
    }
  }else if(tool_name == "slms_3"){
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

  # I'm pretty sure the below check will always pass, but I'll leave it in here just in case
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
  base_path = check_config_and_value("project_base")
  base_path_remote = check_config_and_value("project_base",config_name="default")

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

  # if maf_data provided, no need to loop over seq_type or get from indexed flatfile
  if(!missing(maf_data)){
    muts_region = maf_data %>%
      dplyr::filter(Chromosome == chromosome & Start_Position > qstart & Start_Position < qend) %>%
      dplyr::filter(Tumor_Sample_Barcode %in% these_samples_metadata$sample_id) %>%
      dplyr::filter(t_alt_count >= min_read_support)

    if(streamlined == TRUE){
      muts_region = muts_region %>% dplyr::select(Start_Position, Tumor_Sample_Barcode)
    }else if(basic_columns == TRUE){
      # would error if provided maf_data does not have these columns
      # but if it's truly maf data, it should
      muts_region = muts_region %>% dplyr::select(names(maf_header)[c(1:45)])
    }
  }else if(tool_name == "strelka2"){ # separated out bc only one seq_type, slms_3 case below loops over seq_type
    maf_partial_path = check_config_and_value("results_flatfiles$ssm$all$strelka2")
    #use glue to get the absolute path
    maf_path = glue::glue(maf_partial_path)
    full_maf_path = paste0(base_path, maf_path)
    if(verbose){
      print(maf_partial_path)
      print(full_maf_path)
    }
    full_maf_path_comp = gsub('.{3}$', 'bed',full_maf_path)
    full_maf_path_comp = paste0(full_maf_path_comp, ".gz")
    #check if file is existing or missing
    if(!file.exists(full_maf_path_comp)){
      print(paste("missing:", full_maf_path_comp))
      check_host(verbose=TRUE)
      if(verbose){
        print("using local file")
        print(paste("HERE:",full_maf_path_comp))
      }
      stop(paste("failed to find the file needed for this:",full_maf_path_comp))
    }
    #get tabix command
    tabix_command = paste(tabix_bin, full_maf_path_comp, region)
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
      muts_region = suppressMessages(vroom::vroom(I(muts),
                                     progress = FALSE,
                                     col_types = paste(maf_column_types, collapse = ""),
                                     col_names = maf_columns,
                                     delim = "\t"))
    }
    # filter to only the samples of interest
    muts_region = dplyr::filter(muts_region, Tumor_Sample_Barcode %in% these_samples_metadata$sample_ids)

    if(verbose){
      print('SUCCESS')
    }
  }else if(tool_name == "slms_3"){ # loop over seq_type and unlist at the end
    if(!all(unique(these_samples_metadata$seq_type) %in% c("genome", "capture"))){
      warning("CAUTION! More seq_types than genome and capture found in these_samples_metadata.
       Only genome and capture will be used to get variants.")
      these_samples_metadata = dplyr::filter(these_samples_metadata, seq_type %in% c("genome", "capture"))
    }
    # set up lists
    seq_type_sample_ids = list() # items are vectors of sample_ids
    for(a_seq_type in unique(these_samples_metadata$seq_type)){
        seq_type_sample_ids[[a_seq_type]]=dplyr::filter(these_samples_metadata,seq_type==a_seq_type) %>%
        pull(sample_id)
    }
    seq_type_muts = list() # items are dataframes of mutations
    seq_type_muts_region = list() # items are dataframes of mutations

    for(a_seq_type in names(seq_type_sample_ids)){
      seq_type = a_seq_type # needed for glue
      if(augmented){
        maf_partial_path = check_config_and_value("results_flatfiles$ssm$template$merged$augmented")
      }else{
        message("augmented set to FALSE. Getting variants from merge.")
        maf_partial_path = check_config_and_value("results_flatfiles$ssm$template$merged$deblacklisted")
      }
      #use glue to get the absolute path
      maf_path = glue::glue(maf_partial_path)
      full_maf_path = paste0(base_path, maf_path)
      if(verbose){
        print(maf_partial_path)
        print(full_maf_path)
      }
      full_maf_path_comp = paste0(full_maf_path, ".bgz")
      #check if file is existing or missing
      if(!file.exists(full_maf_path_comp)){
        print(paste("missing:", full_maf_path_comp))
        check_host(verbose=TRUE)
        if(verbose){
          print("using local file")
          print(paste("HERE:",full_maf_path_comp))
        }
        stop(paste("failed to find the file needed for this:",full_maf_path_comp))
      }
      #get tabix command
      tabix_command = paste(tabix_bin, full_maf_path_comp, region, "| cut -f", paste(maf_indexes, collapse = ","))
      if(verbose){
        print(tabix_command)
      }
      #execute the tabix command
      seq_type_muts[[a_seq_type]] = system(tabix_command, intern = TRUE)
      if(verbose){
        print(paste("TYPES:"))
        print(maf_column_types)
        print("NAMES:")
        print(maf_columns)
        print("NUM:")
        print(length(seq_type_muts[[a_seq_type]]))
      }

      if(length(seq_type_muts[[a_seq_type]])==0){
        maf_types_sep = str_split(maf_column_types, pattern = "")[[1]] %>%
            str_replace_all("c", "character") %>%
            str_replace_all("l|i|n", "numeric")
        if(verbose){
          print("adding dummy line")
        }
        seq_type_muts_region[[a_seq_type]] = read.table(textConnection(""), col.names = maf_columns, colClasses = maf_types_sep)
      }else{
      #this is what gets executed with default parameters
      seq_type_muts_region[[a_seq_type]] = suppressMessages(vroom::vroom(I(seq_type_muts[[a_seq_type]]),
                                     progress = FALSE,
                                     col_types = paste(maf_column_types, collapse = ""),
                                     col_names = maf_columns,
                                     delim = "\t"))
      }

      if(verbose){
        print('SUCCESS')
      }
      if(augmented){
        # drop poorly supported reads but only from augmented MAF
        seq_type_muts_region[[a_seq_type]] = dplyr::filter(seq_type_muts_region[[a_seq_type]], t_alt_count >= min_read_support)
      }
      # filter to only the samples of interest
      sample_ids = seq_type_sample_ids[[a_seq_type]]
      seq_type_muts_region[[a_seq_type]] = dplyr::filter(seq_type_muts_region[[a_seq_type]], Tumor_Sample_Barcode %in% sample_ids)
    }
    # combine list into one
    muts_region = do.call("rbind", seq_type_muts_region)
  }
  if(streamlined){
    muts_region = muts_region %>%
      dplyr::select(Start_Position, Tumor_Sample_Barcode)
  }

  return(muts_region)
}
