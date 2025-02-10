#' @title Get CNV and coding SSM combined status
#' 
#' @description For each specified chromosome region (gene name), return status 1 if the copy number (CN) 
#'   state is non-neutral, *i.e.* different from 2, or if the region contains any coding simple somatic mutation (SSM).
#' 
#' @details The user can choose from which regions are intended to return only copy number variation (CNV) status, 
#'   only coding SSM status, or at least the presence of one of them. This behavior is controlled by the arguments 
#'   `genes_and_cn_threshs` (column `cn_thresh`) and `only_cnv`.
#'   
#'   This function internally calls the `get_cn_states`, `get_ssm_by_samples` and `get_coding_ssm_status`functions. 
#'   Therefore, many of its arguments are assigned to these functions. If needed, see the documentation of these 
#'   functions for more information. 
#'   
#'   In the case of returning NA values, this is due to the `get_cn_segments` function not being able to internally 
#'   return any copy number segments from the specified chromosome region.
#' 
#' @param genes_and_cn_threshs A data frame with columns "gene_id" and "cn_thresh". The "gene_id" column stores 
#'   gene symbols (characters) which determine the regions to return CNV and/or coding SSM status. The "cn_thresh" 
#'   column stores integers that mean the maximum or minimum CN states to return status 1 (contains CNV) for 
#'   its respective gene. If this integer is below 2 (neutral CN state for diploids), it is taken as the maximum 
#'   (gene consider as tumor suppressor); if above 2, it is the minimum (oncogene); if equal to 2, do not consider
#'   CNV to return status.
#' @param these_samples_metadata The metadata for samples of interest to be included in the returned matrix.
#'   Can be created with `get_gambl_metadata` function.
#' @param maf_df Optional data frame containing the coding variants for your samples (i.e. output from `get_all_coding_ssm`)
#' @param only_cnv A vector of gene names indicating the genes for which only CNV status should be considered, 
#'   ignoring SSM status. Set this argument to "all" or "none" (default) to apply this behavior to all or none 
#'   of the genes, respectively.
#' @param genome_build Reference genome build. Possible values are "grch37" (default) or "hg38".

#' @param include_hotspots Logical parameter indicating whether hotspots object should also be tabulated. Default is TRUE.
#' @param review_hotspots Logical parameter indicating whether hotspots object should be reviewed to include 
#'   functionally relevant mutations or rare lymphoma-related genes. Default is TRUE.
#' @param seg_data Optionally provide the function with a data frame of segments that will be used instead of the GAMBL flatfiles
#' @param include_silent Set to TRUE if you want Synonymous mutations to also be considered
#' @param adjust_for_ploidy Set to TRUE to scale CN values by the genome-wide average per sample
#' @param this_seq_type Deprecated
#' 
#' @return A data frame with CNV and SSM combined status.
#' 
#' @import dplyr
#' @export
#'
#' @examples
#' 
#' # Get sample metadata including a mix of seq_type
#' all_types_meta = get_gambl_metadata() %>% 
#'             dplyr::filter(pathology == "BL")
#' dplyr::group_by(all_types_meta, seq_type) %>% 
#'      dplyr::summarize(n=dplyr::n())
#' 
#' # For MYC and SYNCRIP, return CNV and SSM combined status; for MIR17HG, 
#' # return only CNV status; for CCND3 return only SSM status
#' genes_and_cn_threshs = data.frame(
#'   gene_id=c("MYC", "MIR17HG", "CCND3", "SYNCRIP"),
#'   cn_thresh=c(3, 3, 2, 1)
#' )
#' 
#' genome_cnv_ssm_status = get_cnv_and_ssm_status(
#'                            genes_and_cn_threshs,
#'                            dplyr::filter(all_types_meta,seq_type=="genome"),
#'                            only_cnv = "MIR17HG")
#' 
#' print(dim(genome_cnv_ssm_status))    
#'                        
#' all_seq_type_status = get_cnv_and_ssm_status(
#'                            genes_and_cn_threshs,
#'                            all_types_meta,
#'                            only_cnv = "MIR17HG")
#' 
#' colSums(genome_cnv_ssm_status)
#' print(dim(all_seq_type_status))   
#' 
get_cnv_and_ssm_status = function(genes_and_cn_threshs,
                                  these_samples_metadata,
                                  maf_df,
                                  seg_data,
                                  only_cnv = "none",
                                  genome_build = "grch37",
                                  include_hotspots = TRUE,
                                  review_hotspots = TRUE,
                                  adjust_for_ploidy=FALSE,
                                  include_silent=FALSE,
                                  this_seq_type){
  
  # check parameters
  stopifnot('`genes_and_cn_threshs` argument is missing.' = !missing(genes_and_cn_threshs))
  if(!missing(this_seq_type)){
    stop("this_seq_type is deprecated. This is now determined from the metadata provided.")
  }
  stopifnot('`genes_and_cn_threshs` argument must be a data frame with columns "gene_id" (characters) and "cn_thresh" (integers).' = {
    k = class(genes_and_cn_threshs) == "data.frame" &
      all( c("gene_id", "cn_thresh") %in% names(genes_and_cn_threshs) )
    if(k){
      is.character(genes_and_cn_threshs$gene_id) &
        is.numeric(genes_and_cn_threshs$cn_thresh) &
        identical(genes_and_cn_threshs$cn_thresh, trunc(genes_and_cn_threshs$cn_thresh))
    }
  })
  
  stopifnot('`genome_build` argument must be "grch37" or "hg38."' = genome_build %in% c("grch37", "hg38"))
  
  
  stopifnot('`only_cnv` argument must be "none", "all", or a subset of `genes_and_cn_threshs$gene_id`' = {
    only_cnv == "none" |
      only_cnv == "all" |
      all(only_cnv %in% genes_and_cn_threshs$gene_id)
  })
  
  # get gene regions
  my_regions = GAMBLR.utils::gene_to_region(gene_symbol = genes_and_cn_threshs$gene_id,
                                            projection = genome_build,
                                            sort_regions = FALSE)
  
  if(length(my_regions) < nrow(genes_and_cn_threshs)){
    genes_and_cn_threshs = dplyr::filter(genes_and_cn_threshs, gene_id %in% names(my_regions))
  }
  
  ### cnv
  thresh_2 = genes_and_cn_threshs$cn_thresh == 2
  genes_and_cn_threshs_non_neutral = genes_and_cn_threshs[!thresh_2,]
  check_cnv = nrow(genes_and_cn_threshs_non_neutral) > 0
  these_samples_metadata = dplyr::filter(these_samples_metadata,
                                         seq_type != "mrna")
  if(check_cnv){
    # get cn states
    regions=my_regions[genes_and_cn_threshs_non_neutral$gene_id]
    
    names(regions)=genes_and_cn_threshs_non_neutral$gene_id
    if(missing(seg_data)){
      cn_matrix = get_cn_states(
        regions = regions,
        strategy="custom_regions",
        these_samples_metadata = these_samples_metadata,
        adjust_for_ploidy = adjust_for_ploidy
      )
    }else{
      cn_matrix = segmented_data_to_cn_matrix(
        regions=regions,
        strategy="custom_regions",
        these_samples_metadata = these_samples_metadata,
        seg_data = seg_data,
        adjust_for_ploidy=adjust_for_ploidy
      )
     
    }

    cn_matrix = cn_matrix[these_samples_metadata$sample_id,, drop=FALSE]
    
    # get cnv status
    cnv_status = mapply(function(cnstate, thresh){
      if(thresh < 2){
        cnstate <= thresh
      }else if(thresh == 2){
        cnstate != thresh
      }else if(thresh > 2){
        cnstate >= thresh
      }
    }, cn_matrix, genes_and_cn_threshs_non_neutral$cn_thresh, USE.NAMES = TRUE, SIMPLIFY = FALSE) %>% 
      as.data.frame %>% 
      {. * 1}
    if("name" %in% colnames(genes_and_cn_threshs)){
      colnames(cnv_status) = genes_and_cn_threshs$name
    }
    rownames(cnv_status) = rownames(cn_matrix)
    cnv_status[is.na(cnv_status)]=0
    # if only CNV statuses are desired, output them
    if(only_cnv == "all"){
      return(cnv_status)
    }
    
    # add cnv status as zero to genes whose cn threshold is 2
    if(any(thresh_2)){
      cnv_status = genes_and_cn_threshs$gene_id[thresh_2] %>% 
        { matrix(0, nrow = nrow(cnv_status), ncol = length(.), dimnames = list(NULL, .)) } %>% 
        cbind(cnv_status)
    }
    print(head(cnv_status))
    cnv_status = dplyr::select(cnv_status, genes_and_cn_threshs$gene_id)
  }
  
  ### ssm
  # genes to get ssm status
  if(only_cnv == "nome"){
    genes_to_check_ssm = genes_and_cn_threshs$gene_id
  }else{
    genes_to_check_ssm = genes_and_cn_threshs$gene_id [ !(genes_and_cn_threshs$gene_id %in% only_cnv) ]
  }
  
  # get maf data
  if(missing(maf_df)){
    my_maf = get_all_coding_ssm(
      these_samples_metadata = these_samples_metadata,
      projection = genome_build,
      include_silent = include_silent
    ) %>%
      dplyr::filter(Hugo_Symbol %in% genes_to_check_ssm)
  }else{
    my_maf = maf_df
  }
  
  
  # get ssm status
  ssm_status = get_coding_ssm_status(
    gene_symbols = genes_to_check_ssm,
    these_samples_metadata = these_samples_metadata,
    maf_data = my_maf,
    genome_build = genome_build,
    min_read_support = min_read_support_ssm,
    include_hotspots = include_hotspots,
    include_silent = FALSE,
    augmented = augmented
  ) 
 
  ssm_status = ssm_status %>% column_to_rownames("sample_id")
  
  # add missing regions to ssm_status as zero statuses
  missing_regions = !(genes_and_cn_threshs$gene_id %in% names(ssm_status))
  if(any(missing_regions)){
    ssm_status = genes_and_cn_threshs$gene_id[missing_regions] %>% 
      { matrix(0, nrow(ssm_status), length(.), dimnames = list(NULL, .)) } %>% 
      cbind(ssm_status) 
  }
  ssm_status = dplyr::select(ssm_status, genes_and_cn_threshs$gene_id)
  ssm_status = ssm_status[these_samples_metadata$sample_id,, drop=FALSE]
  
  # combine cnv and ssm status
  if(check_cnv){
    (cnv_status | ssm_status) * 1
  }else{
    ssm_status
  }
}
