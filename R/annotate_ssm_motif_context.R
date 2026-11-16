#' @title Annotate mutations target motif
#'
#' @description Checks for the presence of mutations at a given motif. A wrapper around
#' [GAMBLR.utils::annotate_ssm_context()] that finds the reference FASTA from the GAMBLR config
#' (repo_base) before falling back to the BSgenome package.
#'
#' @details In positions that reference allele has been mutated, it will capture (motif length - 1) before and (motif length + 1) alleles after the mutated position.
#' Then, it looks for the presence of motif in the captured sequence and check if the mutation has occurred in the indexed position, it will return 
#' SITE and if the the motif is present, but the mutation is not in the indexed position, it will return MOTIF. In other cases, it will return FALSE.
#' NA will be shown if the mutation is an indel mutation.
#'
#' @param maf MAF data frame (required columns: Reference_Allele, Chromosome, Start_Position, End_Position)
#' @param motif The motif sequence (default is WRCY)
#' @param index Position of the mutated allele in the motif
#' @param genome_build The genome build for the variants you are working with (default is grch37)
#' @param fastaPath Can be a path to a FASTA file
#'
#' @param ... Further arguments for [GAMBLR.utils::annotate_ssm_context()] (e.g. `aid_alleles`, `motif_strand`).
#'
#' @return A data frame with a `seq` column and one column per motif ("SITE", "MOTIF", "FALSE" or NA for indels and multi-base substitutions).
#'
#' @rawNamespace import(IRanges, except = c("start", "end", "merge", "shift", "collapse", "union", "slice", "intersect", "setdiff", "desc", "reduce", "trim"))
#' @rawNamespace import(GenomicRanges, except = c("start", "end", "merge", "shift", "union", "intersect", "setdiff", "reduce", "trim"))
#' @import Rsamtools readr dplyr BSgenome
#' @export
#'
#' @examples
#' my_maf <- get_coding_ssm() %>% 
#'   dplyr::filter(Hugo_Symbol=="BCL2") %>%
#'   dplyr::arrange(Chromosome,Start_Position,Tumor_Sample_Barcode) %>%
#'   head()
#' 
#' annotated = annotate_ssm_motif_context(maf = my_maf,
#'                                        motif = "WRCY")
#' 
#' dplyr::select(annotated,1,5,6,11,13,16,seq,WRCY)
#'
annotate_ssm_motif_context <- function(maf,
                                       motif = "WRCY",
                                       index = 3,
                                       genome_build,
                                       fastaPath,
                                       ...
){
    if(missing(genome_build)) genome_build = NULL
    if(is.null(genome_build) && "maf_data" %in% class(maf)){
        genome_build = get_genome_build(maf)
    }
    if(missing(fastaPath)){
        if(is.null(genome_build)){
            stop("no genome_build information provided or present in maf")
        }
        # local reference from the GAMBLR config, if present; otherwise BSgenome
        base <- check_config_and_value("repo_base")
        candidates <- c(paste0(base, "ref/lcr-modules-references-STABLE/genomes/", genome_build,
                               "/genome_fasta/genome.fa"),
                        paste0(base, genome_build, ".fa"))
        fastaPath <- candidates[file.exists(candidates)][1]
    }else if(!file.exists(fastaPath)){
        message("fastaPath ", fastaPath, " not found; using the BSgenome package instead")
        fastaPath <- NA
    }
    if(is.na(fastaPath)){
        GAMBLR.utils::annotate_ssm_context(maf, motif = motif, index = index,
                                           genome_build = genome_build, ...)
    }else{
        GAMBLR.utils::annotate_ssm_context(maf, motif = motif, index = index,
                                           genome_build = genome_build,
                                           fastaPath = fastaPath, ...)
    }
}
