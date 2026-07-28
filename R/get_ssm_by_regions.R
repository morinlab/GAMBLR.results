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
#' @param this_seq_type Deprecated. Inferred from these_samples_metadata.
#' @param these_sample_ids Deprecated. Inferred from these_samples_metadata.
#'
#' @return Returns a data frame of variants in 3 column format or in MAF-like
#'  format (one row per mutation). The MAF-like format also carries a
#'  maf_seq_type column recording which seq_type (genome/capture) each
#'  variant came from, so genome- and capture-derived rows for the same
#'  sample_id can still be separated after the two are merged together.
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

  # When regions_bed or regions_list isn't specified, and only an object is give
  # this makes sure it's handled appropriately based on it's object type
  if(!missing(regions_bed)){
    genome_build = check_get_projection(list(this_bed=regions_bed),
                                        projection,
                                        custom_error = "Please specify a projection that matches the genome build of regions_bed")

    regions = apply(regions_bed, 1, bed2region)
  }else if(!missing(regions_list) & (sum((c("bed_data", "genomic_data", "data.frame") %in% class(regions_list))) >= 1)){
    # bed data was given but not specified with regions_bed, so it is saved in parameter regions_list
    genome_build = check_get_projection(list(this_bed=regions_list),
                                        projection,
                                        custom_error = "Please specify a projection that matches the genome build of regions_bed")

    regions = apply(regions_list, 1, bed2region)
  }else if(!missing(regions_list) & class(regions_list) == "character"){
    regions = regions_list
  }else{
    stop("You must supply either a regions_list vector or regions_bed object")
  }
  if(verbose){
    print(regions)
  }

  genome_build = projection

  if(missing(maf_data) && tool_name == "slms_3"){
    # FAST PATH: one tabix -R (batch region) pull per seq_type against the
    # tabix-indexed merged flatfiles, instead of one tabix subprocess PER
    # REGION (what the lapply(regions, get_ssm_by_region) loop below does --
    # still used as a fallback for strelka2 / maf_data). tabix already
    # restricts rows to those overlapping ANY requested region, so no
    # interval join is needed for the (more common) non-streamlined, full-MAF
    # output -- that's just union + unique(), matching what the old code did
    # with region_name immediately dropped anyway. Only `streamlined` output
    # needs a per-row region label, recovered via cool_overlaps() -- cheap
    # here because by this point the data is already cut down to just the
    # rows tabix matched, not the full merged MAF.

    #patch due to change in merging strategy (kept in sync with get_ssm_by_region)
    maf_header = c(1:104)
    names(maf_header) = c("Hugo_Symbol","Entrez_Gene_Id","Center","NCBI_Build","Chromosome",
                          "Start_Position","End_Position","Strand","Variant_Classification",
                          "Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
                          "dbSNP_RS","dbSNP_Val_Status","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode",
                          "Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2","Tumor_Validation_Allele1",
                          "Tumor_Validation_Allele2","Match_Norm_Validation_Allele1","Match_Norm_Validation_Allele2",
                          "Verification_Status","Validation_Status","Mutation_Status","Sequencing_Phase",
                          "Sequence_Source","Validation_Method","Score","BAM_File","Sequencer","Tumor_Sample_UUID",
                          "Matched_Norm_Sample_UUID","HGVSc","HGVSp","HGVSp_Short","Transcript_ID","Exon_Number",
                          "t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count","n_alt_count","all_effects",
                          "Allele","Gene","Feature","Feature_type","Consequence","cDNA_position","CDS_position",
                          "Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM","DISTANCE",
                          "STRAND_VEP","SYMBOL","SYMBOL_SOURCE","HGNC_ID","BIOTYPE","CANONICAL","CCDS","ENSP",
                          "SWISSPROT","TREMBL","UNIPARC","RefSeq","SIFT","PolyPhen","EXON","INTRON","DOMAINS","AF",
                          "AFR_AF","AMR_AF","ASN_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","CLIN_SIG","SOMATIC",
                          "PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","IMPACT","PICK",
                          "VARIANT_CLASS","TSL","HGVS_OFFSET","PHENO","MINIMISED","gnomAD_ASJ_AF","gnomAD_EAS_AF",
                          "gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF")

    if(basic_columns){
      maf_columns = names(maf_header)[c(1:45)]
      maf_column_types =  "ciccciiccccccclcccclllllllllllllllccccciiiiii"
    }else{
      maf_columns = names(maf_header)
      maf_column_types = "ciccciiccccccclcccclllllllllllllllccccciiiiiiccccccccccccinnccccccccccccccccccclcccccccccnclcncccclnnnnn"
    }
    maf_indexes = maf_header[maf_columns]
    maf_indexes = maf_indexes[order(maf_indexes)]
    maf_columns = names(maf_indexes)
    maf_indexes = unname(maf_indexes)

    # these_samples_metadata was already filtered (seq_type != "mrna",
    # excluded samples dropped) earlier in this function; no need to redo it.
    if(!all(unique(these_samples_metadata$seq_type) %in% c("genome", "capture"))){
      warning("CAUTION! More seq_types than genome and capture found in these_samples_metadata. Only genome and capture will be used.")
      these_samples_metadata = dplyr::filter(these_samples_metadata, seq_type %in% c("genome", "capture"))
    }

    tabix_bin = check_config_and_value("dependencies$tabix")
    base_path = check_config_and_value("project_base")

    # Parse "chr:start-end" -> chrom/start/end once; reused below both for the
    # -R regions file and (if streamlined) for the region-label join.
    regions_split = strsplit(regions, ":")
    pos_split = strsplit(vapply(regions_split, `[`, character(1), 2), "-")
    regions_parsed = data.frame(
      chrom = vapply(regions_split, `[`, character(1), 1),
      start = as.numeric(vapply(pos_split, `[`, character(1), 1)),
      end   = as.numeric(vapply(pos_split, `[`, character(1), 2)),
      stringsAsFactors = FALSE
    )

    # tabix -R expects a TAB-delimited CHROM/POS/POS_TO file (1-based,
    # inclusive) for non-.bed files -- NOT "chr:start-end" per line (confirmed
    # empirically: colon-dash format silently matches zero rows, since with no
    # tabs present each line is read as a single literal chromosome name).
    regions_file = tempfile(fileext = ".txt")
    write.table(regions_parsed, regions_file, sep = "\t", quote = FALSE,
               row.names = FALSE, col.names = FALSE)
    on.exit(unlink(regions_file), add = TRUE)

    seq_type_sample_ids = list()
    for(a_seq_type in unique(these_samples_metadata$seq_type)){
        seq_type_sample_ids[[a_seq_type]] = dplyr::filter(these_samples_metadata, seq_type == a_seq_type) %>%
          pull(sample_id)
    }

    seq_type_muts_region = list()
    for(a_seq_type in names(seq_type_sample_ids)){
      seq_type = a_seq_type # needed for glue
      if(augmented){
        maf_partial_path = check_config_and_value("results_flatfiles$ssm$template$merged$augmented")
      }else{
        maf_partial_path = check_config_and_value("results_flatfiles$ssm$template$merged$deblacklisted")
      }
      maf_path = glue::glue(maf_partial_path)
      full_maf_path_comp = paste0(base_path, maf_path, ".bgz")
      if(!file.exists(full_maf_path_comp)){
        print(paste("missing:", full_maf_path_comp))
        check_host(verbose = TRUE)
        stop(paste("failed to find the file needed for this:", full_maf_path_comp))
      }

      # ONE tabix call for ALL regions (batch/-R mode), instead of one per region
      tabix_command = paste(tabix_bin, "-R", regions_file, full_maf_path_comp,
                            "| cut -f", paste(maf_indexes, collapse = ","))
      if(verbose) print(tabix_command)
      muts = system(tabix_command, intern = TRUE)

      if(length(muts) == 0){
        maf_types_sep = str_split(maf_column_types, pattern = "")[[1]] %>%
          str_replace_all("c", "character") %>% str_replace_all("l|i|n", "numeric")
        seq_type_muts_region[[a_seq_type]] = read.table(textConnection(""), col.names = maf_columns,
                                                        colClasses = maf_types_sep)
      }else{
        seq_type_muts_region[[a_seq_type]] = suppressMessages(vroom::vroom(
          I(muts), progress = FALSE, col_types = paste(maf_column_types, collapse = ""),
          col_names = maf_columns, delim = "\t"))
      }
      if(augmented){
        seq_type_muts_region[[a_seq_type]] = dplyr::filter(seq_type_muts_region[[a_seq_type]],
                                                           t_alt_count >= min_read_support)
      }
      sample_ids = seq_type_sample_ids[[a_seq_type]]
      seq_type_muts_region[[a_seq_type]] = dplyr::filter(seq_type_muts_region[[a_seq_type]],
                                                         Tumor_Sample_Barcode %in% sample_ids)
      # Stamp seq_type onto the rows before the rbind below merges genome and
      # capture pulls together -- same maf_seq_type convention already used
      # by get_ssm_by_sample(), get_ssm_by_genes(), and get_coding_ssm(), so a
      # caller (or a downstream consumer like GAMBLR.data's bundled-data
      # build) can still separate genome-derived from capture-derived rows
      # for a sample_id that has both, instead of the two becoming
      # indistinguishable once merged.
      seq_type_muts_region[[a_seq_type]] = dplyr::mutate(seq_type_muts_region[[a_seq_type]],
                                                         maf_seq_type = a_seq_type)
    }
    muts_all = do.call("rbind", seq_type_muts_region)

    if(streamlined){
      # Need a per-row region label; cheap now since muts_all only contains
      # rows tabix already matched to the requested regions. Reuses the
      # chrom/start/end already parsed above for the -R regions file.
      regions_df = regions_parsed
      regions_df$name = if(use_name_column) regions_bed[["name"]] else regions
      region_mafs = GAMBLR.helpers::cool_overlaps(
        muts_all, regions_df,
        columns2 = c("chrom", "start", "end")
      ) %>%
        dplyr::mutate(start = Start_Position, sample_id = Tumor_Sample_Barcode, region_name = name) %>%
        dplyr::select(start, sample_id, region_name)
    }else{
      region_mafs = muts_all %>%
        unique() %>%
        GAMBLR.utils::create_maf_data(., genome_build = genome_build)
      return(region_mafs)
    }
  }else{
    # fallback (unchanged): strelka2 tool_name, or maf_data provided directly
    # -- neither is the region-count-driven bottleneck this fast path targets
    if(missing(maf_data)){
        region_mafs = lapply(regions, function(x){get_ssm_by_region(
        region = x,
        these_samples_metadata = these_samples_metadata,
        streamlined = streamlined,
        basic_columns = basic_columns,
        tool_name = tool_name,
        augmented = augmented,
        projection = genome_build,
        min_read_support = min_read_support,
        verbose = verbose
        )})
    }else{
      region_mafs = lapply(regions, function(x){get_ssm_by_region(
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
        )})
    }

    if(!use_name_column){
      names(region_mafs) <- regions
    }else{
      names(region_mafs) <- regions_bed[["name"]]
    }

    region_mafs <- list_rbind(region_mafs, names_to = "region_name")
    if(streamlined){
      region_mafs = mutate(region_mafs, start = Start_Position, sample_id =Tumor_Sample_Barcode) %>%
        dplyr::select(start, sample_id, region_name)
    }else{
      region_mafs = region_mafs %>%
        select(-region_name) %>%
        unique() %>%
        GAMBLR.utils::create_maf_data(., genome_build=genome_build)
      return(region_mafs)
    }
  }
}
