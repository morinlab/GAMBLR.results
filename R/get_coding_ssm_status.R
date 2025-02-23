#' @title Get Coding SSM Status.
#'
#' @description Tabulate mutation status (SSM) for a set of genes.
#'
#' @details This function takes a data frame (in MAF-like format) and converts
#' it to a binary one-hot encoded matrix of mutation status for either a set of
#' user-specified genes (via gene_symbols) or, if no genes are provided, default
#' to all lymphoma genes. The default behaviour is to assign each gene/sample_id
#' combination as mutated only if there is a protein coding mutation for that
#' sample in the MAF but this can be configured to use synonymous variants in
#' some (via include_silent_genes) or all (via include_silent) genes.
#' This function also has other filtering and convenience parameters giving
#' the user full control of the return. For more information, refer to the
#' parameter descriptions and examples.
#' Is this function not what you are looking for? Try one of the following,
#' similar, functions; [GAMBLR.results::get_coding_ssm],
#' [GAMBLR.results::get_ssm_by_patients], [GAMBLR.results::get_ssm_by_sample],
#' [GAMBLR.results::get_ssm_by_samples], [GAMBLR.results::get_ssm_by_region],
#' [GAMBLR.results::get_ssm_by_regions].
#'
#' @param gene_symbols A vector of gene symbols for which the mutation status
#'      will be tabulated. If not provided, lymphoma genes will be returned by
#'      default.
#' @param these_samples_metadata The metadata for samples of interest to be
#'      included in the returned matrix. Only the column "sample_id" is
#'      required. If not provided, the matrix is tabulated for all available
#'      samples as default.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original
#'      MAF from each sample for multi-sample patients instead of the augmented
#'      MAF.
#' @param min_read_support Only returns variants with at least this many reads
#'      in t_alt_count (for cleaning up augmented MAFs).
#' @param maf_path If the status of coding SSM should be tabulated from a custom
#'      maf file, provide path to the maf in this argument. The default is set
#'      to NULL.
#' @param maf_data Either a maf loaded from disk or from the database using a
#'      get_ssm function.
#' @param include_hotspots Logical parameter indicating whether hotspots object
#'      should also be tabulated. Default is TRUE.
#' @param keep_multihit_hotspot Logical parameter indicating whether to keep the
#'      gene annotation as mutated when the gene has both hot spot and
#'      non-hotspot mutation. Default is FALSE. If set to TRUE, will report the
#'      number of non-hotspot mutations instead of tabulating for just mutation
#'      presence.
#' @param recurrence_min Integer value indicating minimal recurrence level.
#' @param projection Specify projection (grch37 or hg38) of mutations. Default
#'      is grch37.
#' @param review_hotspots Logical parameter indicating whether hotspots object
#'      should be reviewed to include functionally relevant mutations or rare
#'      lymphoma-related genes. Default is TRUE.
#' @param genes_of_interest A vector of genes for hotspot review. Currently only
#'      FOXO1, MYD88, and CREBBP are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF
#'      file. The default is hg19 genome build.
#' @param include_silent Logical parameter indicating whether to include silent
#'      mutations into coding mutations. Default is FALSE.
#' @param include_silent_genes Optionally, provide a list of genes for which the
#'      Silent variants to be considered. If provided, the Silent variants for
#'      these genes will be included regardless of the include_silent argument.
#' @param suffix Optionally provide a character that will be appended to the end of each name
#' @param this_seq_type Deprecated. This is now determined from the metadata provided.
#'
#' @return A data frame with tabulated mutation status.
#'
#' @import dplyr tidyr GAMBLR.helpers GAMBLR.utils
#' @export
#'
#' @examples
#' # FL Tier 1 and 2 genes
#' genes = dplyr::filter(GAMBLR.data::lymphoma_genes,
#'                       FL==TRUE) %>% dplyr::pull(Gene)
#' 
#' # Metadata for FL genomes and exomes
#' fl_meta = suppressMessages(get_gambl_metadata()) %>% 
#'     dplyr::filter(pathology=="FL",
#'                   cohort != "FL_Crouch",
#'                   seq_type != "mrna")
#' 
#' table(fl_meta$seq_type)
#' # Here, we let the function load the data for us
#'  coding_tabulated_df = get_coding_ssm_status(
#'   gene_symbols=genes,
#'   include_hotspots=FALSE,
#'   genome_build = "hg38",
#'   these_samples_metadata = fl_meta
#'  )
#'  dim(coding_tabulated_df)
#'  head(colnames(coding_tabulated_df))
#' 
#' # Alternatively, we can provide the MAF data directly
#' 
#' # Load the MAF data (let's use the other genome build this time)
#' 
#' maf_data = get_all_coding_ssm(these_samples_metadata = fl_meta,
#'                               projection = "grch37")
#' 
#' coding_tabulated2 = get_coding_ssm_status(gene_symbols=genes,
#'                                           these_samples_metadata = fl_meta,
#'                                           maf_data = maf_data,
#'                                           include_hotspots=FALSE)
#'  dim(coding_tabulated2)
#'  head(colnames(coding_tabulated2))
#' 
get_coding_ssm_status = function(
    gene_symbols,
    these_samples_metadata,
    augmented = TRUE,
    min_read_support = 3,
    maf_path = NULL,
    maf_data,
    include_hotspots = TRUE,
    keep_multihit_hotspot = FALSE,
    recurrence_min = 5,
    review_hotspots = TRUE,
    genes_of_interest = c("FOXO1", "MYD88", "CREBBP"),
    genome_build,
    include_silent = FALSE,
    include_silent_genes,
    suffix,
    this_seq_type
){
  if(!missing(this_seq_type)){
    stop("this_seq_type is deprecated. This is now determined from the metadata provided.")
  }
  if(missing(gene_symbols)){
    message("defaulting to all lymphoma genes")
    gene_symbols = pull(GAMBLR.data::lymphoma_genes, Gene)
  }

  if(!missing(include_silent_genes)){
    message(
        strwrap(
            prefix = " ",
            initial = "",
            "Output will include all genes specified in gene_symbols
            and include_silent_genes parameters."
        )
    )
    gene_symbols <- c(
        gene_symbols,
        include_silent_genes
    ) %>%
    unique()
  }

  if(missing(these_samples_metadata)){
    these_samples_metadata = suppressMessages(get_gambl_metadata())
  }

  if(include_silent){
    message("Including Synonymous variants for all genes...")
    coding_class <- coding_class
  }else{
    coding_class <- coding_class[!grepl("Silent", coding_class)]
  }

  # call it once so the object can be reused if user wants to annotate hotspots
  if(!missing(maf_data)){
    if(!missing(genome_build) & "maf_data" %in% class(maf_data)){
      if(!genome_build == get_genome_build(maf_data)){
        stop("you have specified a genome_build that doesn't match the genome_build attached to maf_data")
      }
    }else if("maf_data" %in% class(maf_data)){
      genome_build = get_genome_build(maf_data)
    }
    coding_ssm = maf_data

  }else if (!is.null(maf_path)){
    message(paste("setting genome_build for MAF to:",genome_build))
    coding_ssm = GAMBLR.utils::create_maf_data(fread_maf(maf_path,genome_build))
  }

  if(missing(maf_data) & is.null(maf_path)){
      coding_ssm = get_all_coding_ssm(projection = genome_build,
                                      these_samples_metadata = these_samples_metadata,
                                      augmented = augmented,
                                      basic_columns = FALSE,
                                      include_silent = include_silent)

  }

  if(missing(include_silent_genes)){
    coding_ssm <- coding_ssm %>%
        dplyr::filter(
            Variant_Classification %in% coding_class
        )
  } else {
    message(
        strwrap(
            prefix = " ",
            initial = "",
            "You have provided gene list with argument include_silent_genes.
            The Silent variants will be included even if the include_silent
            argument is set to FALSE.
            "
        )
    )
    coding_ssm <- coding_ssm %>%
        dplyr::filter(
            Variant_Classification %in% coding_class |
            (
                Hugo_Symbol %in% include_silent_genes &
                Variant_Classification == "Silent"
            )
        )
  }

  coding = coding_ssm %>%
    dplyr::filter(Hugo_Symbol %in% gene_symbols & Variant_Classification != "Synonymous") %>%
    dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    dplyr::rename("sample_id" = "Tumor_Sample_Barcode", "gene" = "Hugo_Symbol") %>%
    unique() %>%
    mutate(mutated = 1)

  samples_table = dplyr::select(these_samples_metadata, sample_id)
  wide_coding = pivot_wider(coding,names_from = "gene", values_from = "mutated",values_fill = 0)
  all_tabulated = left_join(samples_table, wide_coding)
  all_tabulated = all_tabulated %>%
    replace(is.na(.), 0)

  # include hotspots if user chooses to do so
  if(include_hotspots){
    # first annotate
    annotated = annotate_hotspots(coding_ssm, recurrence_min = recurrence_min)
    # review for the supported genes
    if(review_hotspots){
      annotated = review_hotspots(annotated, genes_of_interest = genes_of_interest, genome_build = genome_build)
    }
    message("annotating hotspots")
    hotspots = annotated %>%
      dplyr::filter(Hugo_Symbol %in% genes_of_interest) %>%
      dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, hot_spot) %>%
      dplyr::rename("sample_id" = "Tumor_Sample_Barcode", "gene" = "Hugo_Symbol") %>%
      dplyr::mutate(gene = paste0(gene, "HOTSPOT")) %>%
      unique() %>%
      dplyr::mutate(mutated = ifelse(hot_spot == "TRUE", 1, 0)) %>%
      dplyr::filter(mutated == 1) %>%
      dplyr::select(-hot_spot)

    # long to wide hotspots, samples are tabulated with 0 if no hotspot is detected
    wide_hotspots = pivot_wider(hotspots, names_from = "gene", values_from = "mutated", values_fill = 0)
    # join with the ssm object
    all_tabulated = left_join(all_tabulated, wide_hotspots)
    all_tabulated = all_tabulated %>%
      replace(is.na(.), 0)

    all_tabulated = all_tabulated %>%
      dplyr::select(where(~ any(. != 0)))

    all_tabulated = as.data.frame(all_tabulated)
    # make SSM and hotspots non-redundant by giving priority to hotspot feature and setting SSM to 0
    for (hotspot_site in colnames(wide_hotspots)[grepl("HOTSPOT", colnames(wide_hotspots))]){
      message(hotspot_site)
      this_gene = gsub("HOTSPOT", "", hotspot_site)
      redundant_features = all_tabulated %>%
        dplyr::select(starts_with(this_gene))

      # if not both the gene and the hotspot are present, go to the next iteration
      if(ncol(redundant_features)!= 2) next
      message("OK")
      # if both gene and it's hotspot are in the matrix, give priority to hotspot feature
      all_tabulated[(all_tabulated[, this_gene] >0 & all_tabulated[, paste0(this_gene, "HOTSPOT")] == 1),][,c(this_gene, paste0(this_gene, "HOTSPOT"))][, this_gene] = 0

      # in case gene has both hotspot and another mutation in the same gene,
      # keep both tabulated as multihits
      if(keep_multihit_hotspot){
        # determine which samples have hot spot and another mutation in same gene
        multihits <- annotated %>%
            dplyr::filter(Hugo_Symbol == this_gene) %>%
            group_by(Tumor_Sample_Barcode) %>%
            dplyr::mutate(n_mut = n()) %>%
            dplyr::filter(
                n_mut > 1
            ) %>%
            dplyr::distinct(Tumor_Sample_Barcode, n_mut, hot_spot) %>%
            # account for cases with both hotspot and not hotspot to avoid
            # double-counting the number of mutations
            mutate_at(vars(hot_spot), ~replace_na(., "FALSE")) %>%
            dplyr::mutate(
                n_mut = ifelse(
                    hot_spot == "TRUE",
                    n_mut - 1,
                    n_mut
                )
            ) %>%
            group_by(Tumor_Sample_Barcode) %>%
            dplyr::arrange(n_mut) %>%
            slice_head() %>%
            ungroup %>%
            select(-hot_spot)

        # Return the annotation of this gene to mutated in these samples
        all_tabulated <- all_tabulated %>%
            left_join(
              .,
              multihits,
              by = c("sample_id" = "Tumor_Sample_Barcode")
            ) %>%
            dplyr::mutate(
                {{this_gene}} := ifelse(
                        !is.na(n_mut),
                        n_mut,
                        !!!syms(this_gene)
                    )
            ) %>%
            select(- n_mut)
      }

    }

  }
  if(!missing(suffix)){
    colnames(all_tabulated) = unlist(lapply(colnames(all_tabulated),function(x){paste(x,suffix,sep="_")}))
    colnames(all_tabulated)[1]="sample_id"
  }
  return(all_tabulated)
}
