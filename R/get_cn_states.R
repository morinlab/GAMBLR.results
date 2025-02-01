#' @title Get CN States.
#'
#' @description Get a copy number matrix for all samples based on segmented data in the database.
#'
#' @details This function returns CN states for the specified regions using the CN data from GAMBLR.results and (optionally) assumes regions with no data are diploid.
#' For how to determine/specify the coordinates of each region, refer to the parameter descriptions and examples.
#'
#' @param strategy The general strategy to define regions. Available options are: 'custom_regions','auto_split','cytobands','GISTIC' 
#' @param regions Required when strategy is set to 'custom_regions'. A data frame in bed-like format or a vector of regions in the format "chrom:start-end"
#' @param these_samples_metadata Optional (but highly recommended) metadata table to auto-subset the data to samples in that table before returning. If missing, the result will include a row for every sample in seg_data.
#' @param n_bins_split Split genome into N equally sized bins
#' @param use_cytoband_name Use cytoband names instead of region names, e.g p36.33.
#' @param missing_data_as_diploid Fill in any sample/region combinations with missing data as diploid (e.g., CN state like 2). Default is FALSE.
#' @param projection Specify the genome build you want. Default is grch37.
#' @param adjust_for_ploidy Set to TRUE to scale CN values by the genome-wide average per sample
#'
#' @return Copy number matrix with sample_id as rows and regions as columns.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' 
#' my_metadata = get_gambl_metadata() %>% dplyr::filter(pathology=="BL")
#' 
#' # basic usage with 500 approximately equal-sized bins across the genome
#' # Use the hg38 projection (all coordinates will be relative to hg38)
#' cn_matrix = get_cn_states(strategy="auto_split",
#'                           projection="hg38",
#'                           these_samples_metadata = my_metadata,
#'                           n_bins_split=500)
#'
#' # get the CN state of every lymphoma gene instead
#' 
#' gene_cn = get_cn_states(strategy="custom_regions",
#'                         these_samples_metadata = metadata,
#'                         regions=GAMBLR.data::grch37_lymphoma_genes_bed)
#'                         
#'                         
#' # get the CN state per cytoband using hg38 projection.
#' # adjust CN values to correct for high ploidy
#' cytoband_cn = get_cn_states(strategy="cytobands",
#'                             projection="hg38",
#'                             use_cytoband_name=TRUE,
#'                             these_samples_metadata=metadata,
#'                             adjust_for_ploidy=TRUE)
#'
get_cn_states = function(strategy,
                         regions,
                         #seg_data,
                         these_samples_metadata,
                         n_bins_split = 1000,
                         use_cytoband_name = FALSE,
                         missing_data_as_diploid = FALSE,
                         adjust_for_ploidy=FALSE,
                         projection="grch37"){
  if(!missing(these_samples_metadata)){
    these_samples_metadata = dplyr::filter(these_samples_metadata,seq_type!="mrna")  
    seg_data = get_cn_segments(these_samples_metadata = these_samples_metadata,projection = projection)
  }else{
    seg_data = get_cn_segments(projection = projection)
  }
  
  if(missing(regions)){
    cn_matrix = segmented_data_to_cn_matrix(strategy=strategy,
                                            seg_data=seg_data,
                                            these_samples_metadata = these_samples_metadata,
                                            n_bins_split = n_bins_split,
                                            use_cytoband_name = use_cytoband_name,
                                            missing_data_as_diploid = missing_data_as_diploid,
                                            adjust_for_ploidy = adjust_for_ploidy,
                                            genome_build = projection)
  }else{
    cn_matrix = segmented_data_to_cn_matrix(strategy=strategy,
                                            seg_data=seg_data,
                                            regions=regions,
                                            these_samples_metadata = these_samples_metadata,
                                            n_bins_split = n_bins_split,
                                            use_cytoband_name = use_cytoband_name,
                                            missing_data_as_diploid = missing_data_as_diploid,
                                            adjust_for_ploidy = adjust_for_ploidy,
                                            genome_build = projection)
  }
  return(cn_matrix)
}
