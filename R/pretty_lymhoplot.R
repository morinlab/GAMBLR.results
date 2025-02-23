#' @title LymphGen feature visualization.
#'
#' @description Make a heatmap showing the frequency of LymphGen features across a cohort of samples.
#' 
#' @param these_samples_metadata A data frame containing the metadata for your samples (LymphGen result will be subset to the sample_id in this data frame).
#' @param flavour Required if the output of get_lymphgen has not been provided.
#' @param lymphgen_all The list returned by get_lymphgen.
#' @param with_A53 Set to FALSE if you are using LymphGen results without A53. Recommended option.
#' @param show_side_annotation Set to TRUE if you want the default ComplexHeatmap annotations on the sides of your heatmap. Default is FALSE.
#' @param verbose Set to TRUE to print helpful messages. Useful for debugging. Default is FALSE (not verbose).
#'
#' @return A ComplexHeatmap object
#' @export
#'
#' @import dplyr ggsci stringr tidyr ComplexHeatmap grid GAMBLR.helpers circlize
#'
#' @examples 
#' meta_df = suppressMessages(get_gambl_metadata()) %>% 
#'   dplyr::filter(seq_type %in% c("genome", "capture")) %>%
#'      dplyr::filter(pathology == "DLBCL")
#' pretty_lymphoplot(meta_df, show_side_annotation = TRUE)
#' pretty_lymphoplot(
#'      meta_df,
#'      show_side_annotation = TRUE,
#'      flavour = "with_cnvs.with_sv.no_A53"
#' )
#' 
pretty_lymphoplot = function(
    these_samples_metadata,
    flavour = "with_cnvs.with_sv.with_A53",
    lymphgen_all,
    with_A53 = FALSE,
    show_side_annotation = FALSE,
    verbose = FALSE
){

    if(missing(lymphgen_all)) {
        lymphgen_all <- get_lymphgen(
            flavour = flavour,
            these_samples_metadata = these_samples_metadata,
            keep_original_columns = TRUE
        )
    }

    if(grepl("with_A53",flavour)) {
        with_A53 = TRUE
    }else{  
        message("Running in the mode without A53 ...")
    }

    lymphgen_sample_anno <- lymphgen_all$sample_annotation

    class_ord <- c("EZB", "ST2", "MCD", "BN2", "N1", "A53")

    if(!with_A53){
        class_ord <- class_ord[!class_ord == "A53"]
    }

    full_class_ord <- c(
        class_ord,
        "EZB-COMP", "ST2-COMP", "MCD-COMP", "BN2-COMP", "N1-COMP"
    )

    lymphgen_all$feature_annotation$Class <- factor(
        lymphgen_all$feature_annotation$Class,
        levels = class_ord
    )
    lymphgen_all$lymphgen$LymphGen <- factor(
        lymphgen_all$lymphgen$LymphGen,
        levels = full_class_ord
    )
    lymphgen_feature_anno <- lymphgen_all$feature_annotation %>% 
        dplyr::arrange(Class, desc(prevalence))
    lymphgen_feature_matrix <- lymphgen_all$features %>% 
        column_to_rownames("sample_id")
    lymphgen_classes <- dplyr::select(
        lymphgen_all$lymphgen,
        sample_id,
        LymphGen,
        Model
    ) %>% 
        dplyr::arrange(LymphGen)

    #order the columns (features) based on which class(es) they relate to
    feature_order <- unique(lymphgen_feature_anno$Feature)
    if(verbose){
        print(feature_order)
    }

    #order the rows (samples) based on which class they were assigned to
    sample_ord <- lymphgen_classes$sample_id


    lg_cols <- GAMBLR.helpers::get_gambl_colours("lymphgen")

    lymphgen_feature_matrix <- lymphgen_feature_matrix[
        sample_ord,
        feature_order
    ]
    lymphgen_feature_matrix_char <- lymphgen_feature_matrix
    lymphgen_feature_matrix_char[] <- ""

    #colour each feature per the class it is associated with
    for(feature_row in c(1:nrow(lymphgen_feature_anno))) {
        feature <- lymphgen_feature_anno[feature_row, "Feature"] %>% pull()
        this_class <- lymphgen_feature_anno[feature_row, "Class"] %>% pull()
        for(sample_id in rownames(lymphgen_feature_matrix)) {
            if(
                lymphgen_feature_matrix[sample_id, feature] == 1 &&
                lymphgen_feature_matrix_char[sample_id, feature] == ""
            ){
                lymphgen_feature_matrix_char[
                    sample_id,
                    feature
                ] = as.character(this_class)
            }
        }
    }

    if(show_side_annotation){
        col_fun_EZB <- colorRamp2(c(0, 8), c("white", lg_cols["EZB"]))
        col_fun_ST2 <- colorRamp2(c(0, 8), c("white", lg_cols["ST2"]))
        col_fun_MCD <- colorRamp2(c(0, 10), c("white", lg_cols["MCD"]))
        col_fun_BN2 <- colorRamp2(c(0, 8), c("white", lg_cols["BN2"]))
        col_fun_N1 <- colorRamp2(c(0, 2), c("white", lg_cols["N1"]))
        col_fun_A53 <- colorRamp2(c(0, 8), c("white", lg_cols["A53"]))
        colnames(lymphgen_sample_anno)[1] = "sample_id"

        lymphgen_sample_anno <- left_join(
            lymphgen_classes,
            lymphgen_sample_anno
        )

        #order the rows(samples) based on feature count for each class
        if(with_A53) {
            sample_annotation <- dplyr::arrange(
                lymphgen_sample_anno,
                LymphGen,
                EZB,
                ST2,
                MCD,
                BN2,
                N1,
                A53
            ) %>% column_to_rownames("sample_id")
        } else {
            sample_annotation <- dplyr::arrange(
                lymphgen_sample_anno,
                LymphGen,
                EZB,
                ST2,
                MCD,
                BN2,
                N1
            ) %>% column_to_rownames("sample_id")
        }

        sample_ord <- rownames(sample_annotation)
        sample_annotation <- sample_annotation[
            sample_ord,
            c("LymphGen", class_ord)
        ]

        row_anno <- HeatmapAnnotation(
            df = sample_annotation,
            which = "row",
            col = list(
                LymphGen = lg_cols,
                EZB = col_fun_EZB,
                ST2 = col_fun_ST2,
                MCD = col_fun_MCD,
                BN2 = col_fun_BN2,
                N1 = col_fun_N1,
                A53 = col_fun_A53
            ),
            show_legend = FALSE
        )

        feature_annotation <- group_by(
            lymphgen_feature_anno,
            Feature
        ) %>% 
        slice_head(n = 1) %>%
        ungroup() %>%
        column_to_rownames("Feature") 

        feature_annotation <- feature_annotation[
            feature_order,
            c("prevalence", "proportion_in")
        ]
        colnames(feature_annotation) <- c("incidence", "specificity")
        column_anno <- HeatmapAnnotation(
            df = feature_annotation,
            which = "column",
            col = list(
                Class = lg_cols
            ),
            show_legend = FALSE
        )

        hm <- Heatmap(
            lymphgen_feature_matrix_char[sample_ord,],
            cluster_columns = FALSE,
            cluster_rows = FALSE,
            col = lg_cols,
            column_names_gp = grid::gpar(fontsize = 5),
            left_annotation = row_anno,
            bottom_annotation = column_anno,
            na_col = "white",
            show_row_names = FALSE,
            show_heatmap_legend = FALSE
        )
        draw(hm)
        return(hm)
    }
    hm <- Heatmap(
        lymphgen_feature_matrix_char,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        na_col = "white",
        col = lg_cols,
        column_names_gp = grid::gpar(fontsize = 5),
        show_row_names = FALSE)
    draw(hm)
    return(hm)
}
