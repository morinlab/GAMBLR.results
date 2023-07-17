#' @title Get Lymphgen.
#'
#' @description Get a specific flavour of LymphGen from the main GAMBL outputs.
#'
#' @details Get a specific flavour of LymphGen from the main GAMBL outputs and tidy the composites.
#' Optionally return a matrix of features instead
#'
#' @param flavour Lymphgen flavour.
#' @param these_samples_metadata A metadata table to auto-subset the data to samples in that table before returning.
#' @param return_feature_matrix Boolean parameter, default is FALSE.
#' @param return_feature_annotation Boolean parameter, default is FALSE.
#' @param lymphgen_file Path to lymphgen file.
#' @param keep_all_rows Boolean parameter, default is FALSE.
#' @param keep_original_columns Boolean parameter, default is FALSE.
#'
#' @return A data frame.
#'
#' @import config dplyr tidyr readr stringr tibble
#' @export
#'
#' @examples
#' \dontrun{
#' lymphgens = get_lymphgen(flavour = "no_cnvs.no_sv.with_A53")
#' }
#'
get_lymphgen = function(these_samples_metadata,
                        flavour,
                        return_feature_matrix = FALSE,
                        return_feature_annotation = FALSE,
                        lymphgen_file,
                        keep_all_rows = FALSE,
                        keep_original_columns = FALSE){

  if(missing(these_samples_metadata)){
    if(!keep_all_rows){
      these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome")
    }
  }
  if(missing(flavour)){
    if(!missing(lymphgen_file)){
      lg_path = lymphgen_file
    }else{
      message("please provide a path to your lymphgen output file or one of the following flavours")
      print(check_config_value(config::get("results_merged_wildcards")$lymphgen_template))
      return(NULL)
    }
  }else{
    lg_path = paste0(check_config_value(config::get("project_base")), check_config_value(config::get("results_merged")$lymphgen_template))
    lg_path = glue::glue(lg_path)
  }

  lg = suppressMessages(read_tsv(lg_path))
  lg_tidy = tidy_lymphgen(lg,lymphgen_column_in = "Subtype.Prediction",lymphgen_column_out = "LymphGen")
  if(return_feature_matrix | return_feature_annotation){
    lg_ord = select(lg_tidy,Sample.Name,LymphGen) %>% arrange(LymphGen) %>% pull(Sample.Name)
    lg_levels = select(lg_tidy,Sample.Name,LymphGen) %>% arrange(LymphGen) %>% pull(LymphGen)
    all_mcd = separate(lg_tidy,col="MCD.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>% unique()
    all_mcd_genes = str_remove(all_mcd,"_.*")%>% unique()
    all_mcd_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_mcd_genes)
    feat_mcd = separate(lg_tidy,col="MCD.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1)

    feat_mcd_genes = separate(lg_tidy,col="MCD.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()


    mcd_mat = left_join(all_mcd_df,feat_mcd_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_mcd = mutate(feat_mcd_genes,Class="MCD")

    all_ezb = separate(lg_tidy,col="EZB.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>% pull(MCD) %>% unique()
    all_ezb_genes = str_remove(all_ezb,"_.*")%>% unique()

    feat_ezb_genes = separate(lg_tidy,col="EZB.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()

    all_ezb_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_ezb_genes)
    feat_ezb = separate(lg_tidy,col="EZB.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1)

    ezb_mat = left_join(all_ezb_df,feat_ezb_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_ezb = mutate(feat_ezb_genes,Class="EZB")

    all_bn2 = separate(lg_tidy,col="BN2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>% pull(MCD) %>% unique()
    all_bn2_genes = str_remove(all_bn2,"_.*")%>% unique()
    all_bn2_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_bn2_genes)

    feat_bn2 = separate(lg_tidy,col="BN2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name,Feature) %>% mutate(present=1)

    feat_bn2_genes = separate(lg_tidy,col="BN2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()


    bn2_mat = left_join(all_bn2_df,feat_bn2_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_bn2 = mutate(feat_bn2_genes,Class="BN2")

    all_st2 = separate(lg_tidy,col="ST2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>% unique()
    all_st2_genes = str_remove(all_st2,"_.*") %>% unique()
    all_st2_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_st2_genes)

    feat_st2 = separate(lg_tidy,col="ST2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name,Feature) %>% mutate(present=1)

    feat_st2_genes = separate(lg_tidy,col="ST2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()


    st2_mat = left_join(all_st2_df,feat_st2_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_st2 = mutate(feat_st2_genes,Class="ST2")

    all_n1 = separate(lg_tidy,col="N1.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>% pull(MCD) %>% unique()
    all_n1_genes = str_remove(all_n1,"_.*") %>% unique()

    #all_n1_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_n1)
    all_n1_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_n1_genes)
    feat_n1 = separate(lg_tidy,col="N1.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name,Feature) %>% mutate(present=1)

    feat_n1_genes = separate(lg_tidy,col="N1.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()

    #n1_mat = left_join(all_n1_df,feat_n1) %>% mutate(present=replace_na(present,0)) %>%
    #  pivot_wider(names_from="Feature",values_from="present")

    n1_mat = left_join(all_n1_df,feat_n1_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_n1 = mutate(feat_n1_genes,Class="N1")

    all_genes = c(all_n1_genes,all_ezb_genes,all_st2_genes,all_bn2_genes,all_mcd_genes)
    print(table(all_genes))
    feat_all = bind_rows(feat_n1,feat_st2,feat_mcd,feat_ezb,feat_bn2)

  if(return_feature_annotation){
    #just give the user the association between each feature and its class along with some summary stats
    feat_count = group_by(feat_all,Feature,Class) %>% count()
    return(feat_count)
  }

  all_mat = left_join(ezb_mat,mcd_mat)
  all_mat = left_join(all_mat,bn2_mat)
  all_mat = left_join(all_mat,n1_mat)
  all_mat = left_join(all_mat,st2_mat)
  if(!keep_all_rows){
    all_mat = dplyr::filter(all_mat,Sample.Name %in% these_samples_metadata$sample_id)
  }
  all_mat = all_mat %>% column_to_rownames("Sample.Name")

    return(all_mat)
  }else{
    if(!keep_original_columns){
      lg_tidy = dplyr::select(lg_tidy,Sample.Name,LymphGen)
    }
    lg_tidy = lg_tidy %>% dplyr::rename("sample_id"="Sample.Name")
    if(!keep_all_rows){
      lg_tidy = dplyr::filter(lg_tidy,sample_id %in% these_samples_metadata$sample_id)
    }
    return(lg_tidy)

  }
}
