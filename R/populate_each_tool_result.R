#' @title Populate Each Tool Results.
#'
#' @description Convenience function for returning results from a specified tool.
#'
#' @details This function takes a tool name `tool` as well as other parameters for specifying the requested result.
#' Other parameters include `genome_build`, this can be just one parameter or a vector with different genome builds to return results for.
#' Similarly, `unix_group` can take either one value or a vector with all the different unix groups to return results for.
#' Lastly, the user can subset the returned results to only silent mutations.
#' This is done with setting `include_silent = TRUE` (default is FALSE).
#'
#' @param tool Name of tool to get results from.
#' @param genome_builds A single genome build or a vector of all genome builds to process.
#' @param unix_groups A single unix group or a vector of all unix groups to process.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE.
#'
#' @return Nothing.
#'
#' @import purrr dplyr RMariaDB DBI tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' tool_results = populate_each_tool_result("smlims_3", "grch37", "gambl", FALSE)
#'
#' @keywords internal
populate_each_tool_result = function(tool,
                                     genome_builds,
                                     unix_groups,
                                     include_silent = FALSE){

  database_name = GAMBLR.helpers::check_config_value(config::get("database_name"))
  con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  all_meta = get_gambl_metadata()
  generic_update = function(field_name, sample_id, field_value){
    #note: we'll need to handle strings differently here once we start adding them
    for(i in c(1:length(field_value))){
      fv = field_value[i]
      sid = sample_id[i]
      if(is.numeric(fv)){
        update_q = paste0("UPDATE derived_data set ", field_name, " = ", fv, " WHERE sample_id = \"", sid, "\";")
        print(update_q)
      }else{
        #need to add quotes to the value also
        update_q = paste0("UPDATE derived_data set ", field_name, " = \"", fv, "\" WHERE sample_id = \"", sid, "\";")
        print(update_q)
      }
      dbExecute(con, update_q)
    }
  }
  #check if we're missing sample_ids from sample_table
  #sample_table = config::get("tables")$samples
  #derived_table = config::get("tables")$derived
  #sample_ids = pull(sample_table,sample_id)
  #for(id in sample_ids){
  #  check_q = paste0("select count(*) from ",derived_table, " where sample_id = \"",id,"\";")
  #  num = dbGetQuery(con, check_q)
  #  if(num ==0){
  #    insert_q = paste0("insert into ", derived_table, " (sample_id) values(\"",id,"\");")
  #    print(insert_q)
  #    dbExecute(con, insert_q)
  #  }
  #}
  if(tool == "QC"){
    #flag any cases with QC issues raised
    collated = collate_results()

    qc_issues = dplyr::select(collated,sample_id,QC_flag) %>%
      dplyr::filter(!is.na(QC_flag)) %>%
      dplyr::mutate(flagged = "yes")

    generic_update(sample_id = qc_issues$sample_id, field_name = "QC_issue", field_value = qc_issues$QC_flag)
  }
  if(tool == "sequenza"){
    parse_sequenza = function(sequenza_files){

      seq_data=sequenza_files %>%
        purrr::map(suppressMessages(read_tsv)) %>% #read each file into a list of tibbles
        purrr::map(head, 1) %>% #just keep the first line
        purrr::reduce(rbind) %>% #rbind the elements all back into one
        dplyr::rename(sequenza_cellularity = cellularity, sequenza_ploidy = ploidy) #change the column names

      return(seq_data)
    }

    #seq_files_gambl_hg38 = fetch_output_files(genome_build = "hg38", base_path = "gambl/sequenza_current", results_dir = "02-sequenza", tool_name = "sequenza")
    #sequenza_results = parse_sequenza(seq_files_gambl_hg38$full_path)
    #generic_update(sample_id = seq_files_gambl_hg38$tumour_sample_id, field_name = "sequenza_purity", field_value = sequenza_results$sequenza_cellularity)
    #generic_update(sample_id = seq_files_gambl_hg38$tumour_sample_id, field_name = "sequenza_ploidy", field_value = sequenza_results$sequenza_ploidy)
    #seq_files_gambl = fetch_output_files(genome_build = "grch37", base_path = "gambl/sequenza_current", results_dir = "02-sequenza", tool_name = "sequenza")
    #sequenza_results = parse_sequenza(seq_files_gambl$full_path)
    #generic_update(sample_id = seq_files_gambl$tumour_sample_id, field_name = "sequenza_purity", field_value = sequenza_results$sequenza_cellularity)
    #generic_update(sample_id = seq_files_gambl$tumour_sample_id, field_name = "sequenza_ploidy", field_value = sequenza_results$sequenza_ploidy)

    for(unix_group in unix_groups){
      for(genome_build in genome_builds){
        seq_files_gambl = fetch_output_files(build = genome_build, unix_group = unix_group, base_path = paste0(unix_group, "/sequenza_current"), results_dir = "02-sequenza", tool = "sequenza")
        sequenza_results = parse_sequenza(seq_files_gambl$full_path)
        generic_update(sample_id = seq_files_gambl$tumour_sample_id, field_name = "sequenza_purity", field_value = sequenza_results$sequenza_cellularity)
        generic_update(sample_id = seq_files_gambl$tumour_sample_id, field_name = "sequenza_ploidy", field_value = sequenza_results$sequenza_ploidy)
      }
    }

    #seq_files_gambl_hg38 = fetch_output_files(genome_build = "hg38", base_path = "icgc_dart/sequenza_current", results_dir = "02-sequenza", tool_name = "sequenza")
    #sequenza_results = parse_sequenza(seq_files_gambl_hg38$full_path)
    #generic_update(sample_id = seq_files_gambl_hg38$tumour_sample_id, field_name = "sequenza_ploidy", field_value = sequenza_results$sequenza_ploidy)
    #seq_files_gambl_grch37 = fetch_output_files(genome_build = "hs37d5", base_path = "icgc_dart/sequenza_current", results_dir = "02-sequenza", tool_name = "sequenza")
    #sequenza_results = parse_sequenza(seq_files_gambl_grch37$full_path)
    #generic_update(sample_id = seq_files_gambl_grch37$tumour_sample_id, field_name = "sequenza_ploidy", field_value = sequenza_results$sequenza_ploidy)
  }
  if(tool == "manta"){
    #fetch the output file names per group/build combination
    message("processing results from manta")
    message(paste(unix_groups, sep = ","))
    #separately process by unix group
    for(ug in unix_groups){
      #files_df = find_files_extract_wildcards(tool_name = "manta", genome_build = genome_builds, search_pattern = ".bed", unix_group = ug)
      files_df = find_expected_outputs(tool_name = "manta", unix_group = ug)
      print(head(files_df))
      message(paste("processing", ug))

      n_missing =  files_df %>%
        dplyr::filter(is.na(file_timestamp)) %>%
        count() %>%
        pull(n)

      if(n_missing){
        message(paste("missing outputs for", n_missing))

        files_df=files_df %>%
          dplyr::filter(!is.na(file_timestamp))
      }
      #are there unexpected outputs?
      dupes = names(which(table(files_df$tumour_sample_id) > 1))
      if(length(dupes) > 0){
        message("DUPLICATE RESULTS EXIST. PLEASE FIX THIS AND RERUN")
        message(dupes)
        return()
      }

      all_tsb = pull(all_meta, sample_id)
      no_meta = unique(pull(files_df[which(!files_df$tumour_sample_id %in% all_tsb),], tumour_sample_id))
      if(length(no_meta) > 0){
        message("DROPPING RESULTS FROM SAMPLES WITH NO METADATA:")
        print(no_meta)

        files_df = files_df %>%
          dplyr::filter(!tumour_sample_id %in% no_meta)
      }
      #TODO: Flag and drop files with too many SVs before merging (i.e. remove really bad data). Done manually currently.
      manta_df = process_all_manta_bedpe(files_df, group = ug) #need to add this to the database. Not currently automated

      manta_df %>%
        group_by(tumour_sample_id) %>%
        tally() %>%
        arrange(desc(n)) #have a look at the top offenders (most SVs)
    }
  }
  if(tool == "slms3"){
    gambl_mut_maf = tbl(con, "maf_slms3_hg19_icgc")
    #additional bookkeeping: set matched/unmatched information in the analysis table based on the matched normal ID
    gambl_mutation_normals = gambl_mut_maf %>%
      dplyr::select(Tumor_Sample_Barcode) %>%
      group_by(Tumor_Sample_Barcode) %>%
      as.data.frame()

    gambl_meta_normals = get_gambl_metadata(tissue_status_filter = c('tumour', 'normal')) %>%
      dplyr::select(patient_id, sample_id, tissue_status) %>%
      pivot_wider(id_cols = patient_id, names_from = tissue_status, values_from = sample_id)
    #the above has tumour and normal as separate columns for all paired samples.

    gambl_normals = mutate(gambl_normals, slms3_pairing_status = case_when())
    #just use the mutation table to get summary counts per sample and add to the derived table for convenience

    #update for gambl cases then do the same for icgc, then repeat for coding changes
    gambl_counts = gambl_mut_maf %>%
      dplyr::select(Tumor_Sample_Barcode) %>%
      group_by(Tumor_Sample_Barcode) %>%
      tally() %>%
      as.data.frame()

    generic_update(sample_id = gambl_counts$Tumor_Sample_Barcode, field_name = "slms3_ssm_total", field_value = gambl_counts$n)

    #gambl_vafs = gambl_maf %>%
    #   group_by(Tumor_Sample_Barcode) %>%
    #   mutate(vaf = mean(t_alt_count / (t_ref_count + t_alt_count))) %>%
    #   select(Tumor_Sample_Barcode, vaf) %>%
    #   as.data.frame()

    vaf_q = "select Tumor_Sample_Barcode, avg(t_alt_count/(t_ref_count + t_alt_count)) as vaf from maf_slms3_hg19_icgc group by Tumor_Sample_Barcode"

    vaf_tbl = dbGetQuery(con,vaf_q)

    #update all vafs at once
    generic_update(sample_id = vaf_tbl$Tumor_Sample_Barcode, field_name = "slms3_mean_vaf", field_value = vaf_tbl$vaf)

    if(!include_silent){
      coding_class = coding_class[coding_class != "Silent"]
    }

    #this is where things get SLOW! It seems the most efficient here to do this query per sample.
    for(sam in gambl_counts$Tumor_Sample_Barcode){
      #for(sam in some){
      #check here to see if there's a non-null value and skip if possible
      check_q = paste0("select count(*) as n from derived_data where sample_id = \"", sam, "\" and slms3_ssm_coding is not NULL;")

      num = dbGetQuery(con, check_q) %>%
        pull(n)

      print(paste(sam, num))
      if(num == 0){
        print(paste("working on:", sam))

        coding_num = gambl_mut_maf %>%
          dplyr::filter(Tumor_Sample_Barcode == sam & Variant_Classification %in% coding_class) %>%
          dplyr::count() %>%
          dplyr::pull(n)

        generic_update(sample_id = sam, field_name = "slms3_ssm_coding", field_value = coding_num)
      }else{
        print(paste0("skipping ", sam))
      }
    }
  }
  if(tool == "battenberg_ploidy"){
    # parse purity and ploidy values from copy number caller and add to database
    parse_batt = function(batt_file){
      batt_data =  batt_file %>%
        purrr::map(suppressMessages(read_tsv)) %>%
        purrr::reduce(rbind) %>%
        dplyr::rename(battenberg_cellularity = cellularity, battenberg_ploidy = ploidy, battenberg_psi = psi)

      return(batt_data)
    }

    for(unix_group in unix_groups){
      message(unix_group)
      for(genome_build in genome_builds){
        message(genome_build)
        files = fetch_output_files(build = genome_build, unix_group = unix_group,
                                   base_path = paste0(unix_group, "/battenberg_current"),
                                   results_dir = "02-battenberg", tool = "battenberg_ploidy")

        results_table = files %>%
          mutate(parse_batt(full_path))

        generic_update(sample_id = results_table$tumour_sample_id, field_name = "battenberg_psi", field_value = results_table$battenberg_psi)
        generic_update(sample_id = results_table$tumour_sample_id, field_name = "battenberg_ploidy", field_value = results_table$battenberg_ploidy)
        generic_update(sample_id = results_table$tumour_sample_id, field_name = "battenberg_purity", field_value = results_table$battenberg_cellularity)
        #results_table = files %>%
        #   dplyr::mutate(parse_batt(full_path))

        #sequenza_results = parse_sequenza(seq_files_gambl$full_path)

        #generic_update(sample_id = seq_files_gambl$tumour_sample_id, field_name = "sequenza_purity", field_value = sequenza_results$sequenza_cellularity)
        #generic_update(sample_id = seq_files_gambl$tumour_sample_id, field_name = "sequenza_ploidy", field_value = sequenza_results$sequenza_ploidy)
      }
    }

    #files_gambl_hg38 = fetch_output_files(genome_build = "hg38", base_path = "gambl/battenberg_current", results_dir = "02-battenberg")
    #table column structure is as follows:
    # {tool}_{variable} e.g. battenberg_ploidy and battenberg_purity

    #files_gambl_grch37 = fetch_output_files(genome_build = "grch37", base_path = "gambl/battenberg_current", results_dir = "02-battenberg")
    #results_table = files_gambl_grch37 %>%
    #   mutate(parse_batt(full_path))

    #generic_update(sample_id = results_table$tumour_sample_id, field_name = "battenberg_psi", field_value = results_table$battenberg_psi)
    #generic_update(sample_id = results_table$tumour_sample_id, field_name = "battenberg_ploidy", field_value = results_table$battenberg_ploidy)
    #generic_update(sample_id = results_table$tumour_sample_id, field_name = "battenberg_purity", field_value = results_table$battenberg_cellularity)
  }
}
