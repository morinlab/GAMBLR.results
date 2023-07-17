#' @title Get Manta SVs
#'
#' @description Retrieve Manta SVs and filter.
#'
#' @details Return Manta SVs with additional VCF information to allow for filtering of high-confidence variants.
#' To return SV calls for multiple samples, give `these_sample_ids` a vector of sample IDs, if only one sample is desired,
#' give this parameter one sample ID, as a string (or a vector of characters). The user can also call the `these_samples_metadata`
#' parameter to make use of an already subset metadata table. In this case, the returned calls will be restricted to the sample_ids
#' within that data frame. This function relies on a set of specific internal functions [GAMBLR::id_ease] and [GAMBLR::get_manta_sv_by_samples] (if `from_cache = FALSE`).
#' This function can also restrict the returned calls to any genomic regions specified within `chromosome`, `qstart`, `qend`, 
#' or the complete region specified under `region` (in chr:start-end format), note that chromosome can be either prefixed or not prefixed.
#' Useful filtering parameters are also available, use `min_vaf` to set the minimum tumour VAF for a SV to be returned and `min_score`
#' to set the lowest Manta somatic score for a SV to be returned. `pair_status` can be used to return variants from either matched or unmatched samples.
#' In addition, the user can chose to return all variants, even the ones not passing the filter criteria. To do so, set `pass = FALSE` (default is TRUE).
#' Is it adviseed to run this function with `from_cache = TRUE` (default) to read manta calls from a previous generated merge (cached result).
#' If set to FALSE in combination with `write_to_file = TRUE`, the function will generate new merged manta calls, if the data access restriction allows it.
#' Note, that if `write_to_file` is set to TRUE, the function autodefaults `from_cache = FALSE` to avoid nonsense parameter combinations.
#' Is this function not what you are looking for? Try one of the following, similar, functions;
#' [GAMBLR::get_combined_sv], [GAMBLR::get_manta_sv_by_sample], [GAMBLR::get_manta_sv_by_samples]
#'
#' @param these_sample_ids A vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata A metadata table to auto-subset the data to samples in that table before returning.
#' @param projection The projection genome build. Default is grch37.
#' @param chromosome Optional, the chromosome you are restricting to (can be prefixed or not prefixed).
#' @param qstart Optional, query start coordinate of the range you are restricting to.
#' @param qend Optional, query end coordinate of the range you are restricting to.
#' @param region Optional, region formatted like chrX:1234-5678 (chromosome can be prefixed or not prefixed) instead of specifying chromosome, start and end separately.
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Default is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned. Default is 40.
#' @param pass If TRUE (default) only return SVs that are annotated with PASS in the FILTER column. Set to FALSE to keep all variants, regardless if they PASS the filters.
#' @param pairing_status Use to restrict results (if desired) to matched or unmatched results (default is to return all). This parameter takes the filtering condition as a string ("matched" or "unmatched").
#' @param from_flatfile Set to TRUE by default, FALSE is no longer supported (database).
#' @param verbose Set to FALSE to minimize the output to console. Default is TRUE. This parameter also dictates the verbose-ness of any helper function internally called inside the main function.
#' @param from_cache Boolean variable for using cached results, default is TRUE. If `write_to_file = TRUE`, this parameter auto-defaults to FALSE.
#' @param write_to_file Boolean statement that outputs bedpe file if TRUE, default is FALSE. Setting this to TRUE forces `from_cache = FALSE`.
#'
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs.
#'
#' @import dplyr readr
#' @export
#'
#' @examples
#' #lazily get every SV in the table with default quality filters
#' all_sv = get_manta_sv()
#'
#' #get all SVs for a single sample
#' some_sv = get_manta_sv(these_sample_ids = "94-15772_tumorA")
#'
#' #get the SVs in a region around MYC
#' myc_locus_sv = get_manta_sv(region = "8:128723128-128774067")
#'
get_manta_sv = function(these_sample_ids,
                        these_samples_metadata,
                        projection = "grch37",
                        chromosome,
                        qstart,
                        qend,
                        region,
                        min_vaf = 0.1,
                        min_score = 40,
                        pass = TRUE,
                        pairing_status,
                        from_flatfile = TRUE,
                        verbose = TRUE,
                        from_cache = TRUE,
                        write_to_file = FALSE){
  
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }

  #get samples with the dedicated helper function
  meta_ids = id_ease(these_samples_metadata = these_samples_metadata,
                these_sample_ids = these_sample_ids,
                verbose = verbose,
                this_seq_type = "genome") #only genome samples have manta results

  this_meta = meta_ids$this_metadata
  
  if(write_to_file){
    from_cache = FALSE #override default automatically for nonsense combination of options
  }
  
  if(from_flatfile){
    if(from_cache){
      #get paths and check for file permissions
      output_base = check_config_value(config::get("project_base"))
      output_file = check_config_value(config::get("results_merged")$manta_sv$icgc_dart)
      output_file = paste0(output_base, output_file)
      output_file = glue::glue(output_file)
      
      permissions = file.access(output_file, 4) #check read permissions
      
      if(permissions == -1){
        message("No permission for unix group icgc_dart found, resorting to samples belonging to unix group gambl...")
        output_file = check_config_value(config::get("results_merged")$manta_sv$gambl)
        output_file = paste0(output_base, output_file)
        output_file = glue::glue(output_file)
      }
      if(verbose){
        message(paste0("\nThe cached results were last updated: ", file.info(output_file)$ctime))
        message("\nReading cached results...\n") 
      }
      
      #check for missingness of merged manta results
      if(!file.exists(output_file)){
        print(paste("missing: ", output_file))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
      }
      
      #read merged data
      manta_sv = suppressMessages(read_tsv(output_file)) %>% 
        dplyr::filter(tumour_sample_id %in% this_meta$sample_id, 
                      VAF_tumour >= min_vaf,
                      SCORE >= min_score)
      
      if(verbose){
        no_manta = setdiff(this_meta$sample_id, manta_sv$tumour_sample_id)
        
        if(length(no_manta) > 0){
          message(paste0("No Manta results found for ", length(no_manta), " samples..."))
          print(no_manta)
        } 
      }
      
    }else{
      if(write_to_file){
        #enforce all samples in the altest metadata to be in the merge, if the user decides to overwrite the cached results.
        this_meta = get_gambl_metadata(seq_type_filter = "genome")
      }
      
      #compile the merge based on selected projection (with no filters)
      if(verbose){
        message("\nFrom cache is set to FALSE, this function is now compiling a new merged results file for the selected projection...") 
      }
      
      manta_sv = get_manta_sv_by_samples(these_samples_metadata = this_meta, 
                                         verbose = verbose, 
                                         min_vaf = 0, 
                                         pass = FALSE, 
                                         min_score = 0, 
                                         projection = projection)
      
      #ensure only sample IDs in the full metadata table are kept (i.e if a sample is not in the metadata table, no manta results for any such sample will sneak its way into the merged results file)
      manta_sv = manta_sv %>%
        dplyr::filter(tumour_sample_id %in% this_meta$sample_id)
      
      if(write_to_file){
        #get paths and check for file permissions
        output_base = check_config_value(config::get("project_base"))
        icgc_dart_file = check_config_value(config::get("results_merged")$manta_sv$icgc_dart)
        icgc_dart_file = paste0(output_base, icgc_dart_file)
        icgc_dart_file = glue::glue(icgc_dart_file)
        icgc_dart_folder = gsub(paste0("manta.genome--", projection, ".bedpe"), "", icgc_dart_file)
        
        icgc_permissions = file.access(icgc_dart_folder, 2) #get write permission for the icgc_dart merge (all samples).
        
        if(icgc_permissions == 0){ #get path to gambl samples only merge, if user has acces to the icgc_dart merge.
          gambl_file = check_config_value(config::get("results_merged")$manta_sv$gambl)
          gambl_file = paste0(output_base, gambl_file)
          gambl_file = glue::glue(gambl_file)
          
          #subset icgc_dart to only gambl samples
          gambl_samples = this_meta %>% 
            dplyr::filter(unix_group == "gambl")
          
          gambl_manta_sv = manta_sv %>% 
            dplyr::filter(tumour_sample_id %in% gambl_samples$sample_id)
          
          #write merges to file
          write_tsv(manta_sv, file = icgc_dart_file, append = FALSE)
          write_tsv(gambl_manta_sv, file = gambl_file, append = FALSE)
        }else{
          stop("You do not have the right permissions to write the manta merged files to disk... ")
        }
      }
    }
  }else{
    stop("\nDatabase usage is deprecated, please set from_flatfile to TRUE...")
  }
  
  #deal with chr prefixes based on the selected projection (if return is to be subset to regions...)
  if(!missing(region) || !missing(chromosome)){
    if(projection == "grch37"){
      if(grepl("chr", chromosome)){
        chromosome = gsub("chr", "", chromosome)
      }
    }else if(projection == "hg38"){
      if(!grepl("chr", chromosome)){
          chromosome = paste0("chr", chromosome)
      }
    }
    
    manta_sv = manta_sv %>%
      dplyr::filter((CHROM_A == chromosome & START_A >= qstart & START_A <= qend) | (CHROM_B == chromosome & START_B >= qstart & START_B <= qend))
  }
  
  if(verbose){
    message("\nThe following VCF filters are applied;")
    message(paste0("  Minimum VAF: ", min_vaf))
    message(paste0("  Minimum Score: ", min_score))
    message(paste0("  Only keep variants passing the quality filter: ", pass))
  }
  
  #PASS filter
  if(pass){
    manta_sv = manta_sv %>%
      dplyr::filter(FILTER == "PASS")
  }
  
  #pairing status filter
  if(!missing(pairing_status)){
    if(verbose){
      message(paste0("  Pairing status: ", pairing_status)) 
    }
    manta_sv = manta_sv %>%
      dplyr::filter(pair_status == pairing_status)
  }
  
  #convert to data frame and print some metrics
  manta_sv = as.data.frame(manta_sv)
  
  if(verbose){
    n_variants = nrow(manta_sv)
    unique_samples = unique(manta_sv$tumour_sample_id)
    message(paste0("\nReturning ", n_variants, " variants from ", length(unique_samples), " sample(s)"))
    message("\nDone!") 
  }
  return(manta_sv)
}
