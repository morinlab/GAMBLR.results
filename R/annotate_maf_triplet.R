#' @title Annotate MAF with triplet context
#'
#' @description Give triple sequence of mutated base with its adjacent bases (-1
#' and +1)
#'
#' @details It gives the reference and alternative alleles and looks for the
#' rows of the data frame based on these values for + strand genes and their
#' complement alleles rows for - strand genes, then it can look for the adjacent
#' bases in that mutation position. Also, it can look for all the SNVs in
#' the MAF data frame and provide triple sequences for them (reverse complement
#' sequence for the - strand).
#'
#' @param maf MAF file (required columns: Reference_Allele, Tumor_Seq_Allele2)
#' @param all_SNVs If TRUE, all triplet sequences of single nucleotide 
#'        variants (SNVs) are returned without filtering for specific reference (ref) or 
#'        alternative (alt) alleles. (Default is TRUE). 
#'        When all_SNVs is TRUE, the ref and alt parameters are ignored.
#' @param ref Reference allele (Only relevant when all_SNVs is FALSE; otherwise, this parameter is ignored.)
#' @param alt Alternative allele (Only relevant when all_SNVs is FALSE; otherwise, this parameter is ignored.)
#' @param projection The genome build projection for the variants you are
#'      working with (default is grch37)
#' @param fastaPath The path to the genome FASTA file corresponding 
#'        to the specified genome build. This file is used to extract sequence context 
#'        when no BSgenome package is provided via bsgenome_name. 
#'        - On GSC systems: If not specified, the function attempts to automatically 
#'          infer the FASTA path from the GAMBL configuration.
#'        - On local systems: A valid local path to a FASTA file must be provided.
#' @param bsgenome_name Specifies the name of a BSgenome data package 
#'        to be used for sequence extraction. This parameter overrides both projection 
#'        and fastaPath if provided.
#'        - Format: `BSgenome.<organism>.<provider>.<assembly>[.<masked>]`.
#'        - Example: `"BSgenome.Hsapiens.UCSC.hg38"` for the human UCSC genome build hg38.
#'        - If a masked genome is required, use a name like `"BSgenome.Hsapiens.UCSC.hg38.masked"`.
#' @param pyrimidine_collapse estimates the mutation strand 
#'        using a pyrimidine collapse strategy:
#'        - Reference alleles C or T are interpreted as mutations on the + strand.
#'        - Reference alleles A or G are interpreted as mutations on the - strand.
#'
#' @return A data frame with two to three extra columns, in case
#'      pyrimidine_collapse = FALSE, it will add triple sequence (seq) and the
#'      strand (mutation_strand). When pyrimidine_collapse = T, another column
#'      also will be added to the those extra columns which shows the mutation
#'      (mutation)
#'
#' @import Rsamtools dplyr BSgenome
#' @export
#'
#' @examples
#' maf <- GAMBLR.data::get_coding_ssm(these_sample_id = "DOHH-2")
#'
#' annotate_maf_triplet(maf)
#' annotate_maf_triplet(maf, all_SNVs = FALSE, "C", "T")
#' annotate_maf_triplet(maf, ref = "C", alt = "T", pyrimidine_collapse = TRUE)
#'
#This function gives triple sequence of provided mutated base
annotate_maf_triplet = function(maf,
                                all_SNVs = TRUE,
                                ref,
                                alt,
                                projection = "grch37",
                                fastaPath,
                                bsgenome_name,
                                pyrimidine_collapse = FALSE){
  genome = ""
  bsgenome_loaded = FALSE
  if (projection == "grch37") {
    maf$Chromosome <- gsub("chr", "", maf$Chromosome)
  } else {
    # If there is a mix of prefixed and non-prefixed options
    maf$Chromosome <- gsub("chr", "", maf$Chromosome)
    maf$Chromosome <- paste0("chr", maf$Chromosome)
  }
  # If there is no fastaPath, it will read it from config key
  # Based on the projection the fasta file which will be loaded is different
  if (missing(fastaPath)){
    base <- check_config_value(config::get("repo_base"))
    fastaPath <- paste0(
      base,
      "ref/lcr-modules-references-STABLE/genomes/",
      projection,
      "/genome_fasta/genome.fa"
    )
    if(!file.exists(fastaPath)){
      #try BSgenome
      installed = installed.genomes()
      if(!missing(bsgenome_name)){
        if(bsgenome_name %in% installed){
          genome = getBSgenome(bsgenome_name)

          bsgenome_loaded = TRUE
        }
      }else{
        print(installed)
        stop("Local Fasta file cannot be found. Supply a path to it with fastaPath or use one of the installed BSGenome options using the bsgenome_name parameter")
      }
    }
  }
  # It checks for the presence of a local fastaPath
  if(!bsgenome_loaded){
    # Create a reference to an indexed fasta file.
    if (!file.exists(fastaPath)) {
      stop("Failed to find the fasta file")
    }
    fasta = Rsamtools::FaFile(file = fastaPath)
  }

  # Store the complement of ref and alt alleles
  complement <- c(
    'A'= 'T',
    'T'= 'A',
    'C'= 'G',
    'G'= 'C'
  )
  if(pyrimidine_collapse){
    maf <- maf %>%
        dplyr::mutate(
            mutation_strand=ifelse(
                Reference_Allele %in% c("A","G"),
                    "-",
                        "+"
            )
        ) %>%
        dplyr::mutate(
            mutation = ifelse(
                Reference_Allele %in% c("A","G"),
                     paste0(complement[Reference_Allele],">",complement[Tumor_Seq_Allele2]),
                        paste0(Reference_Allele,">",Tumor_Seq_Allele2)))
  }else{
    maf = dplyr::mutate(
        maf,
        mutation_strand="+")
  }

  CompRef = complement[ref]
  CompAlt = complement[alt]
  # Provide triple sequence for all the SNVs
    # Using BSgenome
  if (all_SNVs == TRUE){
    if(bsgenome_loaded){
      sequences <- maf %>%
        dplyr::mutate(
          seq = ifelse(
            (nchar(maf$Reference_Allele) == 1 &
               nchar(maf$Tumor_Seq_Allele2) == 1
            ),
            as.character(
                Rsamtools::getSeq(
                    genome,
                    maf$Chromosome,
                    start = maf$Start_Position-1,
                    end = maf$Start_Position + 1,
                    strand = mutation_strand
                )
            ),
            "NA"
          )
        )

    }else{
      # Using local Fasta file
      sequences <- maf %>%
        dplyr::mutate(
          seq = ifelse(
            (nchar(maf$Reference_Allele) == 1 &
               nchar(maf$Tumor_Seq_Allele2) == 1
            ),
            as.character(
                Rsamtools::getSeq(
                    fasta,
                    GenomicRanges::GRanges(
                        maf$Chromosome,
                        IRanges(
                            start = maf$Start_Position - 1,
                            end = maf$Start_Position + 1
                        ),
                            strand = maf$mutation_strand
                    )
                )
            ),
            "NA"
          )
        )

    }

  }else{
    # Provide triple sequence of + strand and reverse complement of - strand
    # Mutations on + strand with chosen ref and alt alleles
    # Mutations on - strand with complement ref and alt alleles
    sequences <- maf %>%
      dplyr::mutate(
        seq = ifelse(
          (maf$mutation_strand == "+" &
             maf$Reference_Allele == ref &
             maf$Tumor_Seq_Allele2 == alt
          )|(
            maf$mutation_strand == "-" &
              maf$Reference_Allele == CompRef &
              maf$Tumor_Seq_Allele2 == CompAlt
          ),
          as.character(
              Rsamtools::getSeq(
                 fasta,
                 GenomicRanges::GRanges(
                    maf$Chromosome,
                    IRanges(
                        start = maf$Start_Position - 1,
                        end = maf$Start_Position + 1
                    ),
                        strand = maf$mutation_strand
                 )
             )
          ),
          "NA"
        )
      )

  }
  return(sequences)
}
