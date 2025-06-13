#' @title Get SSM By Regions.
#'
#' @description Efficiently retrieve all mutations across a range of genomic regions.
#'
#' @details This function internally calls [GAMBLR.results::get_ssm_by_region] to retrieve SSM calls for the specified regions.
#' See parameter descriptions for [GAMBLR.results::get_ssm_by_region] for more information on how the different parameters can be called.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR.results::get_coding_ssm],
#' [GAMBLR.results::get_coding_ssm_status], [GAMBLR.results::get_ssm_by_sample],
#' [GAMBLR.results::get_ssm_by_samples], [GAMBLR.results::get_ssm_by_region]
#'
#' @param regions_list Either provide a vector of regions in the chr:start-end format OR
#' @param regions_bed Better yet, provide a bed file with the coordinates you want to retrieve.
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
#' @param use_name_column If TRUE and your bed-format data frame has a name column 
#'  (must be named "name") these can be used to name your regions. To be used with streamlined = TRUE. Default: FALSE.
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
#' @return Returns a data frame of variants in 3 column format or in MAF-like format (one row per mutation).
#'
#' @import tibble dplyr tidyr GAMBLR.utils parallel
#' @export
#'
#' @examples
#' 
#' # Adds column `name` to the bed-format dataframe
#' # by combining "gene" and "region" values sep by "-"
#' regions_bed = GAMBLR.utils::create_bed_data(
#'    GAMBLR.data::grch37_ashm_regions,
#'    fix_names = "concat",
#'    concat_cols = c("gene","region"),sep="-"
#' ) %>% head(20)
#' 
#' DLBCL_meta = suppressMessages(get_gambl_metadata()) %>% 
#'                 dplyr::filter(pathology=="DLBCL", seq_type == "genome")
#' ashm_MAF = get_ssm_by_regions(regions_bed = regions_bed,
#'                              these_samples_metadata = DLBCL_meta,
#'                              streamlined=FALSE)
#' ashm_MAF %>% dplyr::arrange(Start_Position,Tumor_Sample_Barcode) %>%
#'               dplyr::select(Hugo_Symbol,
#'                     Tumor_Sample_Barcode,
#'                     Chromosome,Start_Position,
#'                     Reference_Allele,Tumor_Seq_Allele2)
#'
#'
get_ssm_by_regions = function(regions_list,
                              regions_bed,
                              these_samples_metadata,
                              maf_data,
                              use_name_column = FALSE,
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

  # Also done in get_ssm_by_region, but this allows passing a df 
  # for these_sample_metadata and avoidng a bunch of output messages
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

  if(streamlined){
    message("Streamlined is set to TRUE, this function will disregard anything specified with basic_columns")
    message("To return a MAF with standard 45 columns, set streamlined = FALSE and basic_columns = TRUE")
    message("To return a maf with all (116) columns, set streamlined = FALSE and basic_columns = FALSE")
  }
  bed2region = function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }

  genome_build = NULL
  if(missing(regions_list)){
    if(!missing(regions_bed)){
      genome_build = check_get_projection(list(this_bed=regions_bed),
                                          projection,
                                          custom_error = "Please specify a projection that matches the genome build of regions_bed")
      
      regions = apply(regions_bed, 1, bed2region)
    }else{
      warning("You must supply either regions_list or regions_bed")
    }
  }else{
    regions = regions_list
  }
  if(verbose){
    print(regions)
  }

  if(missing(maf_data)){
    #projection is required if the MAF was not provided
    if(is.null(genome_build)){
      if(missing(projection)){
        stop("projection is required when maf_data is not provided.")
      }else{
        genome_build = projection
      }
    }
    region_mafs = parallel::mclapply(regions, function(x){get_ssm_by_region(
      region = x,          
      these_samples_metadata = these_samples_metadata,
      streamlined = streamlined,
      basic_columns = basic_columns,
      tool_name = tool_name,
      augmented = augmented,
      projection = genome_build,
      min_read_support = min_read_support,
      verbose = verbose
      )},
      mc.cores = 12)
  }else{
    region_mafs = parallel::mclapply(regions, function(x){get_ssm_by_region(
      region = x,
      these_samples_metadata = these_samples_metadata,
      maf_data = maf_data,
      streamlined = streamlined,
      basic_columns = basic_columns,
      tool_name = tool_name,
      augmented = augmented,
      projection = genome_build,
      min_read_support = min_read_support,
      verbose = verbose
      )},
      mc.cores = 12)
  }

  if(!use_name_column){
    rn = regions
  }else{
    rn = regions_bed[["name"]]
  }

  tibbled_data = tibble(region_mafs, region_name = rn)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)
  
  if(streamlined){
    unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
      dplyr::select(start, sample_id, region_name)
      return(unlisted_df)
  }else{
    region_mafs = do.call(bind_rows, region_mafs)
    region_mafs = GAMBLR.utils::create_maf_data(region_mafs, genome_build=genome_build)
    return(region_mafs)
  }

}
