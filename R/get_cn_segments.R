

#' @title Annotate focal and arm-level CNVs
#'
#' @description Convert segmented CN data into arm- and focal-level gain/loss calls using
#' default DLBCL significance peaks (or user-supplied peak files). Produces both raw
#' copy-number matrices and thresholded “GSM” (genome state matrix) scores.
#'
#' @param seg_data Copy-number segments, typically from [get_cn_segments()], with columns
#'   `ID`, `chrom`, `start`, `end`, and `CN` (and optionally `log.ratio`/`adjusted_CN`).
#' @param these_samples_metadata Sample metadata used by [GAMBLR.utils::segmented_data_to_cn_matrix]
#'   to align sample IDs and seq_type.
#' @param focal_peak_file Path to the focal GISTIC peaks file (default
#'   `DLBCL_focal_peaks.18Aug2024.tsv`; hg38 automatically switches to the lifted file).
#' @param broad_peak_file Path to the arm-level GISTIC significance file (default
#'   `DLBCL_broad_significance.18Aug2024.tsv`; hg38 automatically switches to the lifted file).
#' @param single_del_threshold Absolute CN threshold used for a single-copy loss call.
#' @param double_del_threshold Absolute CN threshold used for a homozygous loss call.
#' @param single_amp_threshold Absolute CN threshold used for a single-copy gain call.
#' @param double_amp_threshold Absolute CN threshold used for a high-level gain call.
#' @param adjust_for_ploidy Pass through to [GAMBLR.utils::segmented_data_to_cn_matrix]; set TRUE
#'   to center CN values on each sample’s average ploidy.
#' @param rounded Currently unused (kept for backward compatibility).
#'
#' @return A list with matrices/vectors:
#' \itemize{
#'   \item `gsm`: thresholded focal + arm scores after adjusting focal calls for arm-level background.
#'   \item `not_arm_adjusted`: thresholded scores without focal adjustment.
#'   \item `all_ones`: same as `gsm` but high-level events collapsed to 1 (present/absent).
#'   \item `cnv`: arm-level CN matrix (unthresholded) for significant arms.
#'   \item `focal_cnv`: focal CN matrix (unthresholded) for significant focal peaks.
#'   \item `unthresholded_gsm`: raw CN values for all arm and focal regions.
#'   \item `arm_level_deletions`, `arm_level_gains`, `arm_level_gains_old`: arm-level scores (0/1/2) with
#'     and without length filtering.
#'   \item `focal_deletions`, `focal_gains`: focal scores (0/1/2) before adjustment for arm background.
#' }
#'
#' @details Arm events must span ≥50% of the chromosome arm to count. Focal deletions are
#' adjusted upward by arm losses; focal gains are adjusted downward by arm gains to avoid
#' double-counting. Uses hg38 peak files automatically when `seg_data` reports hg38.
#'
#' @examples
#' \dontrun{
#' segs <- get_cn_segments(these_samples_metadata = metas, projection = "grch37")
#' gsm_out <- annotate_focal_and_arm_level_CNV(segs, metas)
#' head(gsm_out$gsm)
#' }
#' @export 
annotate_focal_and_arm_level_CNV <- function(seg_data,
these_samples_metadata,
  focal_peak_file = "DLBCL_focal_peaks.18Aug2024.tsv",
  broad_peak_file = "DLBCL_broad_significance.18Aug2024.tsv",
  single_del_threshold = 1.3,
  double_del_threshold = 0.9,
  single_amp_threshold = 2.5,
  double_amp_threshold = 3.5,
  adjust_for_ploidy = FALSE,
  rounded = TRUE){
  #minimum CN delta between the segment and the background
  #single_amp_threshold = single_copy_logr_thresh #2*2^single_copy_logr_thresh
  #single_del_threshold = -1 * single_copy_logr_thresh#2*2^-single_copy_logr_thresh
  #double_amp_threshold = double_copy_logr_thresh #2*2^double_copy_logr_thresh
  #double_del_threshold = -1 * double_copy_logr_thresh #2*2^-double_copy_logr_thresh
  #absolute CN thresholds
  #single_amp_threshold = 2*2^single_copy_logr_thresh
  #single_del_threshold = 2*2^-single_copy_logr_thresh
  #double_amp_threshold = 2*2^double_copy_logr_thresh
  #double_del_threshold = 2*2^-double_copy_logr_thresh
  
  
  arm_length_fraction_threshold = 0.5
  seg_data = mutate(seg_data,size = end - start + 1)
  genome_build = get_genome_build(seg_data)
  focals = read_tsv(focal_peak_file)
  
  if(genome_build == "hg38"){
    #use lifted files instead
    broad_peak_file = "DLBCL_broad_significance.hg38.perchr.tsv"
    focal_peak_file = "DLBCL_focal_peaks.hg38.perchr.tsv"
  }
  arm_level_significance = read_tsv(broad_peak_file)
  arm_level_significance = filter(arm_level_significance,significant_deletion + significant_amplification > 0)
  arm_regions = mutate(arm_level_significance,
    chrom = case_when(chrn < 23 ~ as.character(chrn), chrn == as.character(23) ~ "X", TRUE~ "Y")) %>%
    mutate(size = xEnd - start + 1) %>%
    select(-x1:-x2)
  if(genome_build == "hg38"){
    #add missing chr prefix
    arm_regions = mutate(arm_regions,chrom=paste0("chr",chrom))
  }
  arm_regions = mutate(arm_regions,
    region = paste0(chrom,":",as.integer(start+1),"-",as.integer(xEnd))) 
  arm_regions = mutate(arm_regions,arm = str_to_upper(arm))
  
  regions = arm_regions %>% 
    pull(region)

  names(regions)=arm_regions$arm
  message("running segmented_data_to_cn_matrix with custom",length(regions),"custom regions")

  cnv_mat_arm = segmented_data_to_cn_matrix(seg_data = seg_data,
    these_samples_metadata = these_samples_metadata,
    strategy = 'custom_regions', 
    adjust_for_ploidy = adjust_for_ploidy,
    regions=regions,
    rounded = F)

 #count up arm-level events by direction
 
 del_arm = filter(arm_regions, significant_deletion == 1)
 gain_arm = filter(arm_regions, significant_amplification == 1)

 del_arm_mat = cnv_mat_arm[,del_arm$arm]
 del_arm_score = del_arm_mat
 del_arm_score[del_arm_mat>single_del_threshold]=0
 del_arm_score[del_arm_mat<=single_del_threshold ] = 1
 del_arm_score[del_arm_mat<=double_del_threshold ] = 2

 gain_arm_mat = cnv_mat_arm[,gain_arm$arm]
 gain_arm_score = gain_arm_mat

 gain_arm_score[gain_arm_mat < single_amp_threshold]=0
 
 gain_arm_score[gain_arm_mat >= single_amp_threshold ] = 1
 gain_arm_score[gain_arm_mat >= double_amp_threshold] = 2
 
focal_regions = focals %>% 
mutate(chrom = case_when(chr < 23 ~ as.character(chr), chr == as.character(23) ~ "X", TRUE~ "Y"),
    region = paste0(chrom,":",as.integer(peak_start),"-",as.integer(peak_end))) 
  if(genome_build == "hg38"){
    focal_regions = mutate(focal_regions,region = paste0("chr",region))
  }

focal_regions = mutate(focal_regions,Descriptor = gsub("p","P",Descriptor))
focal_regions = mutate(focal_regions,Descriptor = gsub("q","Q",Descriptor))
focal_regions = mutate(focal_regions,Descriptor = gsub(":",".",Descriptor))
#focal_regions = mutate(focal_regions,Descriptor = gsub("q","Q",Descriptor))
  del_focal = filter(focal_regions,type == "DEL")
  gain_focal = filter(focal_regions,type == "AMP")
  regions = focal_regions %>% 
    pull(region)
  #print(focal_regions)
  names(regions)=focal_regions$Descriptor

cnv_mat_focal = segmented_data_to_cn_matrix(seg_data = seg_data,
    these_samples_metadata = these_samples_metadata,
    strategy = 'custom_regions',
    adjust_for_ploidy = adjust_for_ploidy,
    regions=regions,
  rounded = F)

del_focal_mat = cnv_mat_focal[,del_focal$Descriptor]

del_focal_score = del_focal_mat

del_focal_score[del_focal_mat>single_del_threshold]=0
del_focal_score[del_focal_mat<=single_del_threshold ] = 1
del_focal_score[del_focal_mat<=double_del_threshold ] = 2

gain_focal_mat = cnv_mat_focal[,gain_focal$Descriptor]
gain_focal_score = gain_focal_mat

gain_focal_score[gain_focal_mat< single_amp_threshold]=0
gain_focal_score[gain_focal_mat>=single_amp_threshold ] = 1
gain_focal_score[gain_focal_mat>=double_amp_threshold ] = 2
colnames(del_arm_score) = paste0(colnames(del_arm_score),".DEL")
colnames(gain_arm_score) = paste0(colnames(gain_arm_score),".AMP")

#just the arm name for easier matching
colnames(del_arm_mat) = gsub(".DEL","",colnames(del_arm_mat))
colnames(gain_arm_mat) = gsub(".AMP","",colnames(gain_arm_mat))
gain_focal_score_adj = gain_focal_score
gain_focal_mat_adj = gain_focal_mat
del_focal_score_adj = del_focal_score
del_focal_mat_adj = del_focal_mat

for(arm in colnames(del_arm_mat)){
  focal_this_arm= grep(paste0("^",arm),colnames(del_focal_score_adj),value=T)
    if(length(focal_this_arm)>0){
      background = del_arm_mat[[arm]]
      for(event in focal_this_arm){
        foreground = del_focal_mat[[event]]
        adjusted = foreground + background
        adjusted_score= adjusted
        del_focal_mat_adj[[event]] = adjusted
        adjusted_score[adjusted>single_del_threshold]=0
        adjusted_score[adjusted<=single_del_threshold ] = 1
        adjusted_score[adjusted<=double_del_threshold ] = 2
        del_focal_score_adj[[event]] = adjusted_score
        
      }
    }
}
#revise arm-level events based on segment coverage criteria
# Gain/loss segment must cover at least 50% of the arm!
# 

gain_arm_score_orig = gain_arm_score
for(sample in rownames(gain_arm_score)){
  #if(!sample %in% rownames())
  for(gain_arm in colnames(gain_arm_score)){
    if(gain_arm_score[sample,gain_arm]>0){
      #get relevant segments
      this_arm = gsub(".AMP","",gain_arm)
      arm_region = filter(arm_regions,arm==this_arm)
  
      arm_start = pull(arm_region,start)
      arm_end = pull(arm_region,xEnd)
      
      arm_chrom = pull(arm_region,chrom)
      arm_size = pull(arm_region,size)
      segs = filter(seg_data,ID==sample,chrom==arm_chrom)
      segs = filter(segs,(start >= arm_start & end <= arm_end) | #within
                         (start <= arm_start & end >= arm_start)| #encompassing
                         (start >= arm_start & end >= arm_end) | #offset1
                         (start >= arm_start & end <= arm_end)) %>%
                         filter(CN > single_amp_threshold)
      #trim to constrain to this region
      segs = mutate(segs,
        start=ifelse(start< arm_start,arm_start,start),
        end=ifelse(end > arm_end,arm_end,end)
      )
      #print(segs)
      #stop()
      segs = mutate(segs,frac_size = size/arm_size)
      segs_total = sum(segs$frac_size)
      if(any(segs$frac_size > arm_length_fraction_threshold) || segs_total > arm_length_fraction_threshold){
        #qualifies
        #print(segs)
        #print(segs_total)
        #print(any(segs$frac_size > arm_length_fraction_threshold))
      }else{
        print(paste("setting arm-level back to zero",sample,gain_arm))
        print(segs)
        gain_arm_score[sample,gain_arm] = 0
        #print(arm_region)
        #stop()
      }
    }
  }
}
for(arm in colnames(gain_arm_mat)){
  focal_this_arm= grep(paste0("^",arm),colnames(gain_focal_score_adj),value=T)
    if(length(focal_this_arm)>0){
      background = gain_arm_mat[[arm]]
      names(background) = rownames(gain_arm_mat)
      for(event in focal_this_arm){

        foreground = gain_focal_mat[[event]]
        names(foreground) = rownames(gain_arm_mat)
        adjusted = foreground - background
        adjusted_score= adjusted 
        gain_focal_mat_adj[[event]] = adjusted
        adjusted_score[adjusted < single_amp_threshold]=0
        adjusted_score[adjusted >= single_amp_threshold ] = 1
        adjusted_score[adjusted >= double_amp_threshold ] = 2
        gain_focal_score_adj[[event]] = adjusted_score
        if(event == "18Q21.32.AMP"){
          print(foreground["DLBCL11005T"])
          print(background["DLBCL11005T"])
          print(adjusted["DLBCL11005T"])
          print(adjusted_score["DLBCL11005T"])
        }
        
      }
    }
}

gsm = bind_cols(del_arm_score,gain_arm_score,del_focal_score,gain_focal_score) 
gsm_adjusted = bind_cols(del_arm_score,gain_arm_score,del_focal_score_adj,gain_focal_score_adj) 
gsm_1 = gsm_adjusted
gsm_1[gsm_1>1]=1

#colnames(gsm)=paste0("X",colnames(gsm))
return(list(not_arm_adjusted = gsm,
  all_ones = gsm_1,
   gsm = gsm_adjusted,
   cnv=cnv_mat_arm,
   focal_cnv = cnv_mat_focal,
    unthresholded_gsm = bind_cols(del_arm_mat,gain_arm_mat,del_focal_mat,gain_focal_mat),
    arm_level_deletions=del_arm_score,
    arm_level_gains = gain_arm_score,
    arm_level_gains_old = gain_arm_score_orig,
    focal_deletions = del_focal_score,
    focal_gains = gain_focal_score))
}

#' @title Annotate CN segments with cytobands
#'
#' @description Map each copy-number segment to its overlapping cytogenetic band(s)
#' using bundled cytoband tables for grch37 or hg38.
#'
#' @param seg_data Data frame of CN segments with `chrom`, `start`, and `end`
#'   columns (and any other per-segment fields). The genome build is inferred via
#'   `get_genome_build(seg_data)`.
#' @param genome_build Optional build hint (currently unused; kept for compatibility).
#'
#' @return `seg_data` with cytoband annotations from the matching build, including
#'   `cb.chromosome`, `cb.start`, `cb.end`, and `cytoband` (the cytoband name).
#'
#' @details Uses `cytobands_grch37` or `cytobands_hg38` and `cool_overlaps` to join
#' segments to cytobands; overlapping bands are returned as separate rows if applicable.
#'
#' @examples
#' \dontrun{
#' segs <- get_cn_segments(these_samples_metadata = metas)
#' segs_cb <- assign_segment_to_cytoband(segs)
#' head(segs_cb$cytoband)
#' }
#' @export
assign_segment_to_cytoband <- function(seg_data, genome_build){
  gb = get_genome_build(seg_data)
  if(gb=="grch37"){
    arm_df = cytobands_grch37
  }else if(gb == "hg38"){
    arm_df = cytobands_hg38
  }else{
    stop(paste("cannot find arms for",gb))
  }
  matched = cool_overlaps(seg_data,arm_df,
  columns1 = c("chrom","start","end"),
  columns2 =  c("cb.chromosome","cb.start","cb.end"))
  matched = mutate(matched,
    chrom_name = gsub("chr","",chrom),
    cytoband = cb.name
    ) %>% select(-chrom_name,-cb.name)
  matched
}

#' @title Annotate CN segments with chromosome arms
#'
#' @description Map each copy-number segment to its chromosome arm (p/q) for the
#' appropriate genome build (grch37 or hg38).
#'
#' @param seg_data Data frame of CN segments with `chrom`, `start`, and `end`
#'   columns; build is inferred via `get_genome_build(seg_data)`.
#' @param genome_build Optional build hint (currently unused; kept for compatibility).
#'
#' @return `seg_data` with arm annotations, including `chromosome`, `band_start`,
#'   `band_end`, and `arm_name` (e.g., `1p`, `8q`). Overlapping arms (rare) will
#'   yield multiple rows per segment.
#'
#' @details Chooses `chromosome_arms_grch37` or `chromosome_arms_hg38`, then uses
#' `cool_overlaps` to assign segments to arms; `arm_name` is the chromosome plus
#' uppercased arm label.
#'
#' @examples
#' \dontrun{
#' segs <- get_cn_segments(these_samples_metadata = metas)
#' segs_arm <- assign_segment_to_arm(segs)
#' table(segs_arm$arm_name)
#' }
#' @export
assign_segment_to_arm <- function(seg_data,genome_build){
  gb = get_genome_build(seg_data)
  if(gb=="grch37"){
    arm_df = chromosome_arms_grch37 %>% rename(band_start=start,band_end=end)
  }else if(gb == "hg38"){
    arm_df = chromosome_arms_hg38 %>% rename(band_start=start,band_end=end)
  }else{
    stop(paste("cannot find arms for",gb))
  }
  
  matched = cool_overlaps(seg_data,arm_df,
  columns1 = c("chrom","start","end"),
  columns2 =  c("chromosome","band_start","band_end"))
  matched = mutate(matched,
    chrom_name = gsub("chr","",chrom),
    arm_name = paste0(chrom_name,str_to_upper(arm))
    ) %>% select(-chrom_name)
  matched
}

#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL outputs
#'
#' @details This merely loads and returns all the seg_data available for a projection
#' (genome build) and can assign a single value to dummy segments if they are
#' present/identified in the source file
#' @param these_samples_metadata User must provide a metadata table to restrict the data to the samples in your table.
#' The metadata also ensures the proper handling of duplicate sample_id across seq_types and ensures the
#' seq_type in the metadata faithfully represents the seq_type of the data
#' @param flavour Specify what pipeline or source of data to use.
#' Available options are "combined" (for the merge of the best tool for each data type),
#' or one of "purecn_cnvkit", "purecn_denovo", or "battenberg". 
#' Other thabn "combined", the other flavours are incomplete (e.g. limited to genome,
#' matched or capture samples).
#' @param projection Desired genome coordinate system for returned CN segments. Default is "grch37".
#' @param fill_missing_with Specify how to fill values in dummy segments that were created to satisfy GISTIC.
#' The default is "nothing", which causes these to be dropped so empty regions
#' can be handled in subsequent processing steps. For creating a GISTIC input,
#' you would typically want to set this to "avg_ploidy". 
#' This is taken care of for you by [GAMBLR.utils::prepare_gistic_inputs]
#' @param adjust_for_ploidy Set to TRUE to force segments to be adjusted
#' for the average genome-wide copy number (ploidy) of the sample. Beware! There are multiple
#' GAMBLR.* functions that can do this adjustment.
#' If you apply it at this stage you should *NOT* apply it again to the resulting data!
#' 
#' @param verbose Set to TRUE for a chattier experience
#' @param this_seq_type Deprecated.
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr readr glue GAMBLR.utils
#' @export
#'
#' @examples
#' # Example for just exome/capture samples:
#' # Get metadata for just a few capture samples
#' capture_metadata <- suppressMessages(get_gambl_metadata()) %>%
#'   dplyr::filter(seq_type == "capture") %>%
#'   head()
#'
#' # Load the copy number segments for capture samples using hg38 projection
#' capture_segments_hg38 <- get_cn_segments(
#'   these_samples_metadata = capture_metadata,
#'   projection = "hg38"
#' )
#' print(capture_segments_hg38)
#'
#' genome_metadata <- suppressMessages(get_gambl_metadata()) %>%
#'   dplyr::filter(seq_type == "genome") %>%
#'   head()
#' # Create a metadata table with a mix of seq_types
#' mixed_seq_type_meta <- dplyr::bind_rows(capture_metadata, genome_metadata)
#' ## We can load the copy number segments for all samples across seq_types
#' capture_segments_default <- get_cn_segments(
#'   these_samples_metadata = mixed_seq_type_meta
#' )
#' dplyr::group_by(capture_segments_default, ID) %>%
#'   dplyr::summarize(n = dplyr::n())
#' # Note the default projection is "grch37"
#' print(capture_segments_default)
get_cn_segments <- function(these_samples_metadata,
                            projection = "grch37",
                            flavour = "combined",
                            this_seq_type,
                            fill_missing_with = "nothing",
                            adjust_for_ploidy = FALSE,
                            max_CN = 20,
                            verbose = FALSE) {
  if (!missing(this_seq_type)) {
    stop("this_seq_type is deprecated. Subset your metadata instead.")
  }
  if (missing(these_samples_metadata)) {
    message("no metadata provided")
    message("will get segments for all available genome and capture samples")
    these_samples_metadata <- suppressMessages(get_gambl_metadata()) %>% 
      dplyr::filter(seq_type %in% c("genome", "capture"))
    seq_types <- pull(these_samples_metadata, seq_type) %>% unique()
  } else {
    these_samples_metadata <- dplyr::filter(these_samples_metadata,
      seq_type %in% c("genome", "capture"))
    seq_types <- pull(these_samples_metadata, seq_type) %>% unique()
  }

  genome_ids <- dplyr::filter(these_samples_metadata,
    seq_type == "genome") %>% 
    pull(sample_id)
  capture_ids <- dplyr::filter(these_samples_metadata,
    seq_type == "capture") %>% 
    pull(sample_id)
  if (flavour == "combined") {
    cnv_flatfile_template <- check_config_and_value(
      "results_flatfiles$cnv_combined$icgc_dart"
    )
    coltypes = "cciiid"
  } else if (flavour == "battenberg") {
    seq_types = "genome"
    cnv_flatfile_template <- check_config_and_value(
      "results_merged$battenberg"
    )
    coltypes = "cciiidi"
  } else if(grepl("purecn",flavour)){
    seq_types = "capture"
    coltypes = "cciidd"
    cnv_flatfile_template <- check_config_and_value(
      paste0("results_merged$",flavour))
  } else{
    stop("currently available flavours: combined, purecn_denovo, purecn_cnvkit or battenberg")
  }
  df_list <- list()
  for (seq_type in seq_types) {
    cnv_path <- glue::glue(cnv_flatfile_template)
    full_cnv_path <- paste0(check_config_and_value("project_base"), cnv_path)
    # check permissions to ICGC data.
    permissions <- file.access(full_cnv_path, 4)
    if (permissions == -1) {
      message(paste("failed loading from", full_cnv_path[1]))
      message("restricting to non-ICGC data")
      cnv_flatfile_template <- check_config_and_value("results_flatfiles$cnv_combined$gambl")
      cnv_path <- glue::glue(cnv_flatfile_template)
      full_cnv_path <- paste0(check_config_and_value("project_base"), cnv_path)
    }

    # check for missingness.
    if (!file.exists(full_cnv_path)) {
      print(paste("missing: ", full_cnv_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    if (seq_type == "capture") {
      if(verbose){
        print(full_cnv_path)
      }
      
      seg <- suppressMessages(read_tsv(full_cnv_path,
                                      col_types = coltypes,
                                      na = c("NA", "NaN"),
                                      progress = FALSE)) %>%
        dplyr::filter(ID %in% capture_ids)
    } else {
      if(verbose){
        print(full_cnv_path)
      }
      
      seg <- suppressMessages(read_tsv(full_cnv_path,
                                       col_types = coltypes,
                                       na = c("NA", "NaN"),
                                       progress = FALSE)) %>%
        dplyr::filter(ID %in% genome_ids)
    }
    seg <- mutate(seg, seg_seq_type = seq_type)
    df_list[[seq_type]] <- seg
  }
  if (any(unique(df_list[["capture"]]$ID) %in% unique(df_list[["genome"]]$ID))) {
    stop("overlapping IDs found!")
  }

  all_segs <- do.call("bind_rows", df_list)
  if(!"log.ratio" %in% colnames(all_segs)){
    all_segs = rename(all_segs,
      c("log.ratio"="seg.mean"))
  }
  all_segs <- dplyr::mutate(all_segs, CN = ifelse(log.ratio == -2,0,2 * 2^log.ratio))

 

  if (adjust_for_ploidy || "dummy_segment" %in% colnames(all_segs)) {
    if(!"dummy_segment" %in% colnames(all_segs)){
      all_segs$dummy_segment = 0
    }
    if (fill_missing_with == "diploid") {
      if(verbose){
        print("Using diploid")
      }
      
      all_segs <- mutate(all_segs,
        CN = ifelse(dummy_segment == 1, 2, CN),
        log.ratio = ifelse(dummy_segment == 1, 0, log.ratio)
      )
    }
    if (adjust_for_ploidy || fill_missing_with == "avg_ploidy") {
      if(verbose){
        print("Calculating sample-wide Avg_ploidy")
      }
      
      real_segs <- dplyr::filter(
        all_segs,
        dummy_segment == 0
      )
      real_segs <- real_segs %>%
        mutate(
          length = end - start + 1,
          CN_seg = CN * length,
          logr_seg = log.ratio * length
        ) %>%
        group_by(ID) %>%
        summarise(
          mean = mean(CN),
          real_mean = sum(CN_seg) / sum(length),
          real_mean_logr = sum(logr_seg) / sum(length)
        ) # actual average per base
      all_segs <- left_join(all_segs, select(real_segs, ID, real_mean, real_mean_logr), by = "ID")
      all_segs = mutate(all_segs,CN=ifelse(CN>max_CN,max_CN,CN))
      all_segs <- mutate(all_segs,
        CN = ifelse(dummy_segment == 1, real_mean, CN),
        log.ratio = ifelse(dummy_segment == 1, real_mean_logr, log.ratio),
        adjusted_CN = 2 + CN - real_mean,
        rounded_CN = round(adjusted_CN, 0)
      )
    }
    #if (fill_missing_with == "nothing") {
    #  #drop dummy segments entirely
    #  all_segs <- dplyr::filter(all_segs, dummy_segment == 0)
    #} else {
    #  stop("fill_missing_with must be 'nothing', 'diploid', or 'avg_ploidy'")
    #}
  }
  #else{
  #  message("dummy segments are not annotated in the inputs")
  #  message("fill_missing_with parameter will be ignored")
  #}
  # return S3 class with CN segments and genome_build
  all_segs <- create_seg_data(all_segs, projection)
  return(all_segs)
}
