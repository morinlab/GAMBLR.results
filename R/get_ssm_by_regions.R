#' @title Get SSM By Regions.
#'
#' @description Efficiently retrieve all mutations across a range of genomic regions.
#'
#' @details This function internally calls [GAMBLR::get_ssm_by_region] to retrieve SSM calls for the specified regions.
#' See parameter descriptions for [GAMBLR::get_ssm_by_region] for more information on how the different parameters can be called.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR::get_coding_ssm],
#' [GAMBLR::get_coding_ssm_status], [GAMBLR::get_ssm_by_patients], [GAMBLR::get_ssm_by_sample],
#' [GAMBLR::get_ssm_by_samples], [GAMBLR::get_ssm_by_region]
#'
#' @param regions_list Either provide a vector of regions in the chr:start-end format OR.
#' @param regions_bed Better yet, provide a bed file with the coordinates you want to retrieve.
#' @param streamlined Return a basic rather than full MAF format, default is TRUE.
#' @param maf_data Use an already loaded MAF data frame.
#' @param use_name_column If your bed-format data frame has a name column (must be named "name") these can be used to name your regions.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details).
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param seq_type The seq_type you want back, default is genome.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Boolean parameter set to FALSE per default. Set to TRUE to return fewer columns.
#'
#' @return Returns a data frame of variants in MAF-like format.
#'
#' @import tibble dplyr tidyr
#' @export
#'
#' @examples
#' #basic usage, adding custom names from bundled ashm data frame
#' regions_bed = dplyr::mutate(grch37_ashm_regions, name = paste(gene, region, sep = "_"))
#'
#' ashm_basic_details = get_ssm_by_regions(regions_bed = regions_bed)
#'
#' full_details_maf = get_ssm_by_regions(regions_bed = regions_bed,
#'                                       basic_columns=T)
#'
get_ssm_by_regions = function(regions_list,
                              regions_bed,
                              streamlined = TRUE,
                              maf_data = maf_data,
                              use_name_column = FALSE,
                              from_indexed_flatfile = TRUE,
                              mode = "slms-3",
                              augmented = TRUE,
                              seq_type = "genome",
                              projection = "grch37",
                              min_read_support = 4,
                              basic_columns = FALSE){


  bed2region = function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }

  if(missing(regions_list)){
    if(!missing(regions_bed)){
      regions = apply(regions_bed, 1, bed2region)
    }else{
      warning("You must supply either regions_list or regions_df")
    }
  }
  if(missing(maf_data)){
    print(regions)
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x,
                                                                streamlined = streamlined,
                                                                from_indexed_flatfile = from_indexed_flatfile,
                                                                mode = mode,
                                                                augmented = augmented,
                                                                seq_type = seq_type,
                                                                projection = projection,
                                                                basic_columns=basic_columns)})
  }else{
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x,
                                                                streamlined = streamlined,
                                                                maf_data = maf_data,
                                                                from_indexed_flatfile = from_indexed_flatfile,
                                                                mode = mode,
                                                                basic_columns=basic_columns)})
  }
  if(!use_name_column){
    rn = regions
  }else{
    rn = regions_bed[["name"]]
  }

  if(basic_columns){
    #this must always force the output to be the standard set.
    #hence, return everything after binding into one data frame
    print("bind_rows")
    return(bind_rows(region_mafs))
  }
  tibbled_data = tibble(region_mafs, region_name = rn)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)
  if(streamlined){
    unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
      dplyr::select(start, sample_id, region_name)

  }else{
    unlisted_df = mutate(unnested_df, Chromosome = region_mafs$Chromosome, End_Position = region_mafs$End_Position, Start_Position = region_mafs$Start_Position, Tumor_Sample_Barcode = region_mafs$Tumor_Sample_Barcode) %>%
      dplyr::select(Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode, region_name)
  }
  return(unlisted_df)
}
