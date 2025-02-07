#' @title Get ASHM Count Matrix.
#'
#' @description Prepare a matrix with one row per sample and one column per
#' region using a set of hypermutated regions.
#'
#' @details Values are the number of mutations in that patient in the region.
#'
#' @param regions_bed A bed file with one row for each region. 
#' The first three columns in this file MUST contain the Chromosome, Start and End position
#' @param maf_data Optionally provide a data frame in the MAF format, otherwise
#'      either GAMBLR.data or GAMBLR.results will be used.
#' @param these_samples_metadata This is used to complete your matrix. All GAMBL
#'      samples of the specified seq_type will be used by default. Provide a data frame with at least
#'      sample_id for all samples if you are using non-GAMBL data.
#' @param this_seq_type The seq_type to return results for. Must be a single value. Only used if no
#'      metadata is provided with these_samples_metadata.
#' @param projection The genome build we are working with
#'
#' @return A data frame with a row for every sample in these_samples_metadata and a column for every region in regions_bed
#'
#' @import dplyr tibble
#' @export
#'
#' @examples
#' 
#'\dontrun{
#'   DLBCL_genome_meta = get_gambl_metadata() %>% dplyr::filter(pathology=="DLBCL")
#' regions_bed <- dplyr::mutate(
#'      GAMBLR.data::grch37_ashm_regions,
#'      name = paste(gene, region, sep = "_")
#' )
#'
#' matrix <- get_ashm_count_matrix(
#'      regions_bed = regions_bed,
#'      this_seq_type = "genome",
#'      these_samples_metadata = DLBCL_genome_meta
#' )
#'}
get_ashm_count_matrix = function(
        regions_bed,
        maf_data,
        these_samples_metadata,
        this_seq_type = "genome",
        projection
    ){
    if(missing(these_samples_metadata)){
        these_samples_metadata <- get_gambl_metadata() %>%
            dplyr::filter(seq_type == this_seq_type) 
    }else{
        these_samples_metadata <- these_samples_metadata %>%
            dplyr::filter(seq_type == this_seq_type) 
    }
    
    if(missing(regions_bed)){
      if(!missing(projection)){
        message(
          paste("Using aSHM regions in",projection,"projection as regions_bed")
        )
        if(projection=="grch37"){
          regions_bed <- create_bed_data(GAMBLR.data::grch37_ashm_regions,
                                         fix_names = "concat",
                                         concat_cols = c("gene","region"),
                                         sep="-")
            
        }else if(projection == "hg38"){
          regions_bed <- create_bed_data(GAMBLR.data::hg38_ashm_regions,
                                         fix_names = "concat",
                                         concat_cols = c("gene","region"),
                                         sep="-")
          
        }else{
          stop("unsupported projection")
        }
      }else{
        stop("either projection or a regions_bed containing a genome_build is required")
      }
    }else{
      if("bed_data" %in% class(regions_bed)){
        if(missing(projection)){
          projection = get_genome_build(regions_bed)
        }else{
          if(!projection == get_genome_build(regions_bed) ){
            stop("genome_build in regions_bed does not match projection!")
          } 
        }
      }
    }
    
    ashm_maf <- get_ssm_by_regions(
        regions_bed = regions_bed,
        streamlined = TRUE,
        maf_data = maf_data,
        use_name_column = TRUE,
        these_samples_metadata=these_samples_metadata,
        projection=projection,
        this_seq_type=this_seq_type
    )

    ashm_counted <- ashm_maf %>%
        group_by(sample_id, region_name) %>%
        tally()

    
    #fill out all combinations so we can get the cases with zero mutations
    eg <- expand_grid(
        sample_id = pull(these_samples_metadata, sample_id),
        region_name = unique(ashm_counted$region_name)
    )
    all_counts <- left_join(eg, ashm_counted) %>%
        mutate(n = replace_na(n, 0)) %>%
        unique() #not sure where the duplicates are coming from but its annoying

    all_counts_wide <- pivot_wider(
        all_counts,
        id_cols = sample_id,
        names_from = region_name,
        values_from = n
    ) %>%
        column_to_rownames(var = "sample_id")

    return(all_counts_wide)
}
