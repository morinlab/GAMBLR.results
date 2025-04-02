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
#' @param regions_list Either provide a vector of regions in the chr:start-end format OR.
#' @param regions_bed Better yet, provide a bed file with the coordinates you want to retrieve.
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' @param streamlined If TRUE (default), only 3 columns will be kept in the maf (start, sample_id and region name). 
#' To return more columns, set this parameter to FALSE, see `basic_column` for more info. 
#' Note, if this parameter is TRUE, the function will disregard anything specified with `basic_columns`.
#' @param maf_data Use an already loaded MAF data frame.
#' @param use_name_column If your bed-format data frame has a name column (must be named "name") these can be used to name your regions.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details).
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param this_seq_type The seq_type you want back, default is genome.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Parameter to be used when streamlined is FALSE. 
#' Set this parameter to TRUE for returning a maf with standard 45 columns, set to FALSE to keep all 116 maf columns in the returned object. 
#' To return all 116 maf columns, set this parameter to FALSE.
#' @param verbose Boolean parameter set to FALSE per default.
#'
#' @return Returns a data frame of variants in MAF-like format.
#'
#' @import tibble dplyr tidyr GAMBLR.utils
#' @export
#'
#' @examples
#' 
#' regions_bed = GAMBLR.utils::create_bed_data(
#'    GAMBLR.data::grch37_ashm_regions,
#'    fix_names = "concat",
#'    concat_cols = c("gene","region"),sep="-"
#' ) %>% head(20)
#' 
#' DLBCL_meta = suppressMessages(get_gambl_metadata()) %>% 
#'                 dplyr::filter(pathology=="DLBCL")
#' ashm_MAF = get_ssm_by_regions(regions_bed = regions_bed,
#'                              these_samples_metadata = DLBCL_meta,
#'                              streamlined=FALSE)
#' ashm_MAF %>% dplyr::arrange(Start_Position,Tumor_Sample_Barcode) %>%
#'               dplyr::select(Hugo_Symbol,
#'                     Tumor_Sample_Barcode,
#'                     Chromosome,Start_Position,
#'                     Reference_Allele,
#'                     Tumor_Seq_Allele2)
#'
#'
get_ssm_by_regions = function(regions_list,
                              regions_bed,
                              these_sample_ids = NULL,
                              these_samples_metadata = NULL,
                              streamlined = TRUE,
                              maf_data = maf_data,
                              use_name_column = FALSE,
                              from_indexed_flatfile = TRUE,
                              mode = "slms-3",
                              augmented = TRUE,
                              this_seq_type = "genome",
                              projection = "grch37",
                              min_read_support = 4,
                              basic_columns = FALSE,
                              verbose = FALSE){
  genome_build = NULL
  if(missing(these_samples_metadata)){
    stop("these_samples_metadata is required")
  }
  #ensure we only process the seq_type specified
  these_samples_metadata = dplyr::filter(these_samples_metadata,
    seq_type == this_seq_type)
  if(streamlined){
    message("Streamlined is set to TRUE, this function will disregard anything specified with basic_columns")
    message("To return a MAF with standard 45 columns, set streamlioned = FALSE and basic_columns = TRUE")
    message("To return a maf with all (116) columns, set streamlined = FALSE and basic_columns = FALSE")
  }
  bed2region = function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }

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
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x,
          
                                                                these_samples_metadata = these_samples_metadata,
                                                                streamlined = streamlined,
                                                                from_indexed_flatfile = from_indexed_flatfile,
                                                                mode = mode,
                                                                augmented = augmented,
                                                                this_seq_type = this_seq_type,
                                                                projection = genome_build,
                                                                basic_columns = basic_columns)})
  }else{
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x,
                    
                                                                these_samples_metadata = these_samples_metadata,
                                                                this_seq_type = this_seq_type,
                                                                streamlined = streamlined,
                                                                maf_data = maf_data,
                                                                from_indexed_flatfile = from_indexed_flatfile,
                                                                mode = mode,
                                                                basic_columns = basic_columns)})
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
