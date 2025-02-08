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
#' @param all_SNVs To give us all the triplet sequences of SNVs and not
#'      specifying any specific ref and alt alleles (default is TRUE)
#' @param ref Reference allele
#' @param alt Alternative allele
#' @param genome_build The genome build for the variants you are
#'      working with (default is to infer it from the MAF)
#' @param fastaPath Can be a path to a FASTA file on a disk. When on GSC,
#'      this is first attempted to be inferred from the gambl reference through
#'      path specified in config. Local files are also accepted as value here.
#' @param pyrimidine_collapse Estimate mutation_strand and
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
#' maf <- GAMBLR.open::get_coding_ssm(projection="hg38") %>% head(n=500)
#' # peek at the data
#' dplyr::select(maf,1:12) %>% head()
#'
#' maf_anno <- annotate_maf_triplet(maf)
#' dplyr::select(maf_anno,1:12,seq) %>% head()
#' # Each mutation is now associated with it's sequence context in the
#' # reference genome in a column named seq
#' \dontrun{
#'   annotate_maf_triplet(maf, all_SNVs = FALSE, "C", "T")
#'   annotate_maf_triplet(maf, ref = "C", alt = "T", pyrimidine_collapse = TRUE)
#' }

annotate_maf_triplet = function(maf,
                                all_SNVs = TRUE,
                                ref,
                                alt,
                                genome_build,
                                fastaPath,
                                pyrimidine_collapse = FALSE){
  genome = ""
  # Get projection from NCBI_Build column of the maf
  if ("NCBI_Build" %in% colnames(maf)) {
    genome_build = tolower(maf$NCBI_Build[1])
  } # If it is not in the maf, look for it in bsgenome_name
  else if (!missing(bsgenome_name)) {
    bsgenome_parts = unlist(strsplit(bsgenome_name, "\\."))
    genome_build = bsgenome_parts[4]
  } # Cannot find it ask for mutating it
  else{
    stop("Please add the 'NCBI_Build' column to your MAF file to specify the genome build.")
  }
  bsgenome_loaded = FALSE

  # If there is no fastaPath, it will read it from config key

  # Based on the genome_build the fasta file which will be loaded is different
  if (missing(fastaPath)){
    if("maf_data" %in% class(maf)){
      genome_build = get_genome_build(maf)
    }
    if(missing(genome_build)){
      stop("no genome_build information provided or present in maf")
    }
    base <- check_config_value(config::get("repo_base"))
    fastaPath <- paste0(
      base,
      "ref/lcr-modules-references-STABLE/genomes/",
      genome_build,
      "/genome_fasta/genome.fa"
    )
    if(!file.exists(fastaPath)){
      #try BSgenome
      installed = installed.genomes()
      if(genome_build=="hg38"){
        bsgenome_name = "BSgenome.Hsapiens.UCSC.hg38"
      }else if(genome_build == "grch37"){
        bsgenome_name = "BSgenome.Hsapiens.NCBI.GRCh37"
      }else{
        stop(paste("unsupported genome:",genome_build))
      }
      if(bsgenome_name %in% installed){
          genome = getBSgenome(bsgenome_name)
          bsgenome_loaded = TRUE
      }else{
        print(installed)
        print(paste("Local Fasta file cannot be found and missing genome_build",bsgenome_name,"Supply a fastaPath for a local fasta file or install the missing BSGenome package and re-run"))
      }
    }
  }
  # It checks for the presence of a local fastaPath
  if(!bsgenome_loaded){
    # Create a reference to an indexed fasta file.
    if (!file.exists(fastaPath)) {
      stop("Failed to find the fasta file and no compatible BSgenome found")
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
