

#' @title Get SSM By Gene
#'
#' @description Retrieve all SSMs from GAMBL for one or more genes
#'
#' @details 
#'
#' @param genes Gene symbol for one or more genes
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#'
#' @rawNamespace import(vroom, except = c("col_skip", "fwf_positions", "default_locale", "date_names_lang", "cols_only", "output_column", "col_character", "col_guess", "spec", "as.col_spec", "fwf_cols", "cols", "col_date", "col_datetime", "locale", "col_time", "cols_condense", "col_logical", "col_number", "col_integer", "col_factor", "fwf_widths", "date_names_langs", "problems", "date_names", "col_double", "fwf_empty"))
#' @import dplyr DBI stringr glue GAMBLR.helpers
#'
#' @examples
#' dlbcl_meta = get_gambl_metadata() %>% 
#'  dplyr::filter(pathology=="DLBCL", seq_type!= "mrna")
#' genes_maf = get_ssm_by_genes(genes = c("EZH2","KMT2D"),these_samples_metadata = dlbcl_meta)
#' genes_maf %>% 
#'  dplyr::group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>% 
#'  count()
#' 
#' @export
get_ssm_by_genes = function(genes,
                           these_samples_metadata = NULL,
                           projection = "grch37",
                           ashm_regions,
                           drop_silent_outside_ashm_regions = FALSE,
                           expand_by = 0) {
    
    all_ssms = list()
    if(!missing(ashm_regions)){
      message("ashm data frame provided, will use the contents instead of bundled data")
      ashm_coords = ashm_regions
      colnames(ashm_coords)[c(1:3)]=c("chrom","ashm_region_start","ashm_region_end")
      if(projection == "grch37"){
        ashm_coords = mutate(ashm_coords,chrom = gsub("chr","",chrom))
      }
    }else if(projection == "grch37"){
      ashm_coords = GAMBLR.data::grch37_ashm_regions %>% 
        mutate(chrom = gsub("chr","",chr_name)) %>%
        rename(c("ashm_region_start"="hg19_start","ashm_region_end"="hg19_end")) %>% select(-chr_name)
    }else{
      ashm_coords = GAMBLR.data::hg38_ashm_regions %>% 
        rename(c("chrom"="chr_name","ashm_region_start"="hg38_start","ashm_region_end"="hg38_end"))
    }
    for(gene in genes){
        #get gene region first
        gene_region = suppressMessages(gene_to_region(gene,projection=projection))
        these_samples_metadata = filter(these_samples_metadata,seq_type %in% c("genome","capture"))
        
        for(s_type in unique(these_samples_metadata$seq_type)){
            seq_type_ssms = get_ssm_by_region(region=gene_region,
                                            basic_columns=TRUE,
                                            streamlined=FALSE,
                                            these_samples_metadata =  filter(these_samples_metadata,seq_type == s_type),
                                            projection=projection) %>%
                                            filter(Hugo_Symbol==gene)
            all_ssms[[paste0(gene,"-",s_type)]] = seq_type_ssms %>% 
              mutate(maf_seq_type = s_type)
            
        }
        
        all_ssm_maf = do.call("bind_rows",all_ssms)
    }
    if(drop_silent_outside_ashm_regions){
      
      ashm_coords = mutate(ashm_coords,
            ashm_region_start = ashm_region_start - expand_by,
            ashm_region_end = ashm_region_end + expand_by
        )
        coding_ssm_maf = filter(all_ssm_maf,Variant_Classification %in% vc_nonSynonymous)
        silent_ssm_maf = filter(all_ssm_maf,!Variant_Classification %in% vc_nonSynonymous)
        silent_ssm_maf = cool_overlaps(silent_ssm_maf,
        ashm_coords,columns1=c("Chromosome","Start_Position","End_Position"),
        columns2=c("chrom","ashm_region_start","ashm_region_end"))
        all_ssm_maf = bind_rows(silent_ssm_maf,coding_ssm_maf)
    }
  return(all_ssm_maf)
}
