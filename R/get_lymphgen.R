#' @title Get Lymphgen.
#'
#' @description Get a specific flavour of LymphGen from the main GAMBL outputs.
#'
#' @details Get a specific flavour of LymphGen from the main GAMBL outputs and tidy the composites.
#' Optionally return a matrix of features instead
#'
#' @param flavour Lymphgen flavour.
#' @param these_samples_metadata A metadata table to auto-subset the data to samples in that table before returning.
#' @param lymphgen_file Path to lymphgen file.
#' @param keep_all_rows Boolean parameter, default is FALSE.
#' @param keep_original_columns Boolean parameter, default is FALSE.
#' @param streamlined Boolean, set to true to get just a data frame with one column for sample_id and one for LymphGen class
#' @param verbose Boolean, set to TRUE to print informational messages. Useful for debugging. Default is FALSE
#'
#' @return If run with A list of data frames with the following names:
#'  lymphgen (a data frame containing the tidy LymphGen output),
#'  features (a binary matrix indicating which patients had each feature),
#'  feature_annotation (a data frame with one row per LymphGen feature reduced to gene or arm, for arm-level events and summary statistics for the feature across the cohort),
#'  features_long (a data frame with one row per LymphGen feature/patient event),
#'  sample_annotation (a data frame with one row per sample and columns indicating the number of features for each LymphGen class in that sample)

#'
#' @import config dplyr tidyr readr stringr tibble
#' @export
#'
#' @examples
#' my_meta <- get_gambl_metadata()
#' lymphgen_all <- get_lymphgen(
#'   flavour = "no_cnvs.no_sv.with_A53",
#'   these = my_meta,
#'   keep_original_columns = TRUE
#' )
#' head(lymphgen_all$features[, c(1:14)])
#'
#' head(lymphgen_all$lymphgen[, c(1:14)])
#'
get_lymphgen <- function(these_samples_metadata,
                         flavour,
                         lymphgen_file,
                         keep_all_rows = FALSE,
                         keep_original_columns = FALSE,
                         streamlined = FALSE,
                         verbose = FALSE) {
  with_A53 <- FALSE
  if (missing(these_samples_metadata)) {
    if (!keep_all_rows) {
      these_samples_metadata <- get_gambl_metadata(seq_type_filter = "genome")
    }
  }
  if (missing(flavour)) {
    if (!missing(lymphgen_file)) {
      if (grepl("with_A53", lymphgen_file)) {
        with_A53 <- TRUE
      } else {
        message("NO A53")
      }
      lg_path <- lymphgen_file
    } else {
      message(paste("please provide a path to your lymphgen output",
        "file or one of the following flavours"))
      opts <- check_config_and_value(
        "results_merged_wildcards$lymphgen_template"
        )
      print(opts)
      return(NULL)
    }
  } else {
    if (grepl("with_A53", flavour) && grepl("with_cnvs", flavour)) {
      with_A53 <- TRUE
    } else {
      message("NO A53")
    }
    lg_path <- paste0(
      check_config_and_value("repo_base"),
      check_config_and_value("results_versioned$lymphgen_template$default")
    )
    lg_path <- glue::glue(lg_path)
    if (verbose) {
      message(paste("loading from:", lg_path))
    }
  }

  lg <- suppressMessages(read_tsv(lg_path,
    progress = FALSE
  ))
  lg_tidy <- tidy_lymphgen(lg,
      lymphgen_column_in = "Subtype.Prediction",
      lymphgen_column_out = "LymphGen")

  if (!keep_all_rows) {
    lg_tidy <- dplyr::filter(lg_tidy,
                            `Sample.Name` %in% these_samples_metadata$sample_id)
  }
  if (!streamlined) {
    lg_ord <- select(lg_tidy, Sample.Name, LymphGen) %>%
      arrange(LymphGen) %>%
      pull(Sample.Name)
    lg_levels <- select(lg_tidy, Sample.Name, LymphGen) %>%
      arrange(LymphGen) %>%
      pull(LymphGen)
    all_mcd <- suppressWarnings(separate(lg_tidy, col = "MCD.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "MCD") %>%
      dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>%
      unique()

    all_mcd_genes <- str_remove(all_mcd, "_.*") %>% unique()
    all_mcd_df <- expand.grid(Sample.Name = unique(lg_tidy$Sample.Name), Feature = all_mcd_genes)
    feat_mcd <- suppressWarnings(separate(lg_tidy, col = "MCD.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1)

    feat_mcd_genes <- suppressWarnings(separate(lg_tidy, col = "MCD.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1) %>%
      mutate(Feature = str_remove(Feature, "_.*")) %>%
      group_by(Sample.Name, Feature) %>%
      slice_head()


    mcd_mat <- left_join(all_mcd_df,
                         feat_mcd_genes,
                         by = c("Sample.Name", "Feature")) %>%
      mutate(present = replace_na(present, 0)) %>%
      pivot_wider(names_from = "Feature", values_from = "present")
    feat_mcd <- mutate(feat_mcd_genes, Class = "MCD")

    all_ezb <- suppressWarnings(
      separate(lg_tidy,
              col = "EZB.Features",
              into = c(paste0("Feature_MCD_",
                              seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "MCD") %>%
      dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>%
      unique()
    all_ezb_genes <- str_remove(all_ezb, "_.*") %>%
      str_remove("Trans.*") %>%
      unique()

    feat_ezb_genes <- suppressWarnings(
      separate(lg_tidy,
              col = "EZB.Features",
              into = c(paste0("Feature_MCD_",
                              seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1) %>%
      mutate(Feature = str_remove(Feature, "_.*")) %>%
      mutate(Feature = str_remove(Feature, "Trans.*")) %>%
      group_by(Sample.Name, Feature) %>%
      slice_head()

    all_ezb_df <- expand.grid(Sample.Name = unique(lg_tidy$Sample.Name), Feature = all_ezb_genes)
    feat_ezb <- suppressWarnings(separate(lg_tidy, col = "EZB.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1)

    ezb_mat <- left_join(all_ezb_df, feat_ezb_genes, by = c("Sample.Name", "Feature")) %>%
      mutate(present = replace_na(present, 0)) %>%
      pivot_wider(names_from = "Feature", values_from = "present")
    feat_ezb <- mutate(feat_ezb_genes, Class = "EZB")

    all_bn2 <- suppressWarnings(separate(lg_tidy, col = "BN2.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "MCD") %>%
      dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>%
      unique()
    all_bn2_genes <- str_remove(all_bn2, "_.*") %>%
      str_remove("Fusion.*") %>%
      unique()
    all_bn2_df <- expand.grid(Sample.Name = unique(lg_tidy$Sample.Name), Feature = all_bn2_genes)

    feat_bn2 <- suppressWarnings(separate(lg_tidy, col = "BN2.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1)

    feat_bn2_genes <- suppressWarnings(separate(lg_tidy, col = "BN2.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1) %>%
      mutate(Feature = str_remove(Feature, "_.*")) %>%
      mutate(Feature = str_remove(Feature, "Fusion.*")) %>%
      group_by(Sample.Name, Feature) %>%
      slice_head()


    bn2_mat <- left_join(all_bn2_df, feat_bn2_genes, by = c("Sample.Name", "Feature")) %>%
      mutate(present = replace_na(present, 0)) %>%
      pivot_wider(names_from = "Feature", values_from = "present")
    feat_bn2 <- mutate(feat_bn2_genes, Class = "BN2")

    all_st2 <- suppressWarnings(separate(lg_tidy, col = "ST2.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "MCD") %>%
      dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>%
      unique()
    all_st2_genes <- str_remove(all_st2, "_.*") %>% unique()
    all_st2_df <- expand.grid(Sample.Name = unique(lg_tidy$Sample.Name), Feature = all_st2_genes)

    feat_st2 <- suppressWarnings(separate(lg_tidy, col = "ST2.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1)

    feat_st2_genes <- suppressWarnings(separate(lg_tidy, col = "ST2.Features", into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1) %>%
      mutate(Feature = str_remove(Feature, "_.*")) %>%
      group_by(Sample.Name, Feature) %>%
      slice_head()


    st2_mat <- left_join(all_st2_df, feat_st2_genes,
                        by = c("Sample.Name", "Feature")) %>%
      mutate(present = replace_na(present, 0)) %>%
      pivot_wider(names_from = "Feature", values_from = "present")
    feat_st2 <- mutate(feat_st2_genes, Class = "ST2")

    all_n1 <- suppressWarnings(separate(lg_tidy,
                              col = "N1.Features",
                              into = c(paste0("Feature_MCD_",
                              seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "MCD") %>%
      dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>%
      unique()
    all_n1_genes <- str_remove(all_n1, "_.*") %>% unique()

    all_n1_df <- expand.grid(Sample.Name = unique(lg_tidy$Sample.Name),
                            Feature = all_n1_genes)
    feat_n1 <- suppressWarnings(separate(lg_tidy,
                                col = "N1.Features",
                                into = c(paste0("Feature_MCD_",
                                seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"),
                   values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1)

    feat_n1_genes <- suppressWarnings(separate(lg_tidy,
                                      col = "N1.Features",
                                      into = c(paste0("Feature_MCD_",
                                      seq(1:18))), sep = ",")) %>%
      pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name, Feature) %>%
      mutate(present = 1) %>%
      mutate(Feature = str_remove(Feature, "_.*")) %>%
      group_by(Sample.Name, Feature) %>%
      slice_head()

    n1_mat <- left_join(all_n1_df, feat_n1_genes,
                        by = c("Sample.Name", "Feature")) %>%
      mutate(present = replace_na(present, 0)) %>%
      pivot_wider(names_from = "Feature", values_from = "present")
    feat_n1 <- mutate(feat_n1_genes, Class = "N1")
    if (with_A53) {
      all_a53 <- suppressWarnings(
        separate(lg_tidy, col = "A53.Features",
        into = c(paste0("Feature_MCD_", seq(1:18))), sep = ",")) %>%
        pivot_longer(starts_with("Feature_"), values_to = "MCD") %>%
        dplyr::filter(!is.na(MCD)) %>%
        pull(MCD) %>%
        unique()
      all_a53_genes <- str_remove(all_a53, "_.*") %>% unique()

      all_a53_df <- expand.grid(Sample.Name = unique(lg_tidy$Sample.Name),
                                Feature = all_a53_genes)
      feat_a53 <- suppressWarnings(
        separate(lg_tidy,
                col = "A53.Features",
                into = c(paste0("Feature_MCD_",
                          seq(1:18))), sep = ",")) %>%
        pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
        dplyr::filter(!is.na(Feature)) %>%
        select(Sample.Name, Feature) %>%
        mutate(present = 1)

      feat_a53_genes <- suppressWarnings(
        separate(lg_tidy,
          col = "A53.Features",
          into = c(paste0("Feature_MCD_",
            seq(1:18))), sep = ",")) %>%
        pivot_longer(starts_with("Feature_"), values_to = "Feature") %>%
        dplyr::filter(!is.na(Feature)) %>%
        select(Sample.Name, Feature) %>%
        mutate(present = 1) %>%
        mutate(Feature = str_remove(Feature, "_.*")) %>%
        group_by(Sample.Name, Feature) %>%
        slice_head()

      a53_mat <- left_join(all_a53_df,
                          feat_a53_genes, by = c("Sample.Name", "Feature")) %>%
        mutate(present = replace_na(present, 0)) %>%
        pivot_wider(names_from = "Feature", values_from = "present")
      feat_a53 <- mutate(feat_a53_genes, Class = "A53")

      all_genes <- c(all_a53_genes,
                     all_n1_genes,
                     all_ezb_genes,
                     all_st2_genes,
                     all_bn2_genes,
                     all_mcd_genes)
      feat_all <- bind_rows(feat_a53,
                            feat_n1,
                            feat_st2,
                            feat_mcd,
                            feat_ezb,
                            feat_bn2)
      all_mat <- left_join(ezb_mat, mcd_mat, by = "Sample.Name")
      all_mat <- left_join(all_mat, bn2_mat)
      all_mat <- left_join(all_mat, n1_mat)
      all_mat <- left_join(all_mat, st2_mat)
      all_mat <- left_join(all_mat, a53_mat)
    } else {
      all_genes <- c(all_n1_genes,
                     all_ezb_genes,
                     all_st2_genes,
                     all_bn2_genes,
                     all_mcd_genes)
      feat_all <- bind_rows(feat_n1, feat_st2, feat_mcd, feat_ezb, feat_bn2)
      all_mat <- left_join(ezb_mat, mcd_mat)
      all_mat <- left_join(all_mat, bn2_mat)
      all_mat <- left_join(all_mat, n1_mat)
      all_mat <- left_join(all_mat, st2_mat)
    }


    feat_count <- feat_all %>%
      unique() %>%
      mutate(total = length(unique(feat_all$Sample.Name))) %>%
      group_by(Feature, Class, total) %>%
      reframe(affected = sum(present), prevalence = 100 * sum(present) / total) %>%
      unique()

    feat_all_expand <- expand.grid(Sample.Name = unique(lg_tidy$Sample.Name), Feature = "X", present = 0, Class = unique(feat_all$Class))
    feat_all_expand <- bind_rows(feat_all_expand, feat_all)

    sample_count <- group_by(feat_all_expand, Class, `Sample.Name`) %>%
      summarise(n = sum(present))
    sample_count_wide <- pivot_wider(sample_count,
                                    id_cols = `Sample.Name`,
                                    names_from = "Class", values_from = "n")


    if (!keep_all_rows) {
      all_mat <- dplyr::filter(all_mat,
                              Sample.Name %in% these_samples_metadata$sample_id)
    }
    if (!keep_original_columns) {
      lg_tidy <- dplyr::select(lg_tidy, Sample.Name, Model.Used, LymphGen)
    }
    lg_tidy <- lg_tidy %>% dplyr::rename("sample_id" = "Sample.Name",
                                         "Model" = "Model.Used")
    if (!keep_all_rows) {
      lg_tidy <- dplyr::filter(lg_tidy,
                               sample_id %in% these_samples_metadata$sample_id)
    }
    # convert everything we return to use sample_id instead of Sample.Name (Caution: this may break some code that relies on this function)
    sample_count_wide <- dplyr::rename(sample_count_wide, c("sample_id" = "Sample.Name"))
    feat_all <- dplyr::rename(feat_all, c("sample_id" = "Sample.Name"))
    all_mat <- dplyr::rename(all_mat, c("sample_id" = "Sample.Name"))


    # drop composites and other for calculating feature enrichment
    lymphgen_nocomp <- dplyr::filter(lg_tidy, LymphGen %in% c("MCD", "EZB", "BN2", "N1", "A53", "ST2")) %>%
      dplyr::select(sample_id, LymphGen)
    features_with_classified <- left_join(lymphgen_nocomp, feat_all, by = "sample_id")

    features_counted <- group_by(features_with_classified, Feature) %>%
      summarise(in_class = sum(Class == LymphGen), out_class = sum(Class != LymphGen)) %>%
      mutate(proportion_in = (in_class + 1) / (in_class + out_class + 1))

    features_annotated <- left_join(feat_count, features_counted, by = c("Feature"))

    # features_annotated= ungroup(features_counted) %>%

    return(list(
      lymphgen = lg_tidy,
      full_lymphgen = lg,
      features = all_mat,
      feature_annotation = features_annotated,
      features_long = feat_all,
      sample_annotation = sample_count_wide
    ))
  } else { # streamlined output
    if (!keep_original_columns) {
      lg_tidy <- dplyr::select(lg_tidy, Sample.Name, LymphGen)
    }
    lg_tidy <- lg_tidy %>% dplyr::rename("sample_id" = "Sample.Name")
    if (!keep_all_rows) {
      lg_tidy <- dplyr::filter(lg_tidy, sample_id %in% these_samples_metadata$sample_id)
    }
    return(lg_tidy)
  }
}
