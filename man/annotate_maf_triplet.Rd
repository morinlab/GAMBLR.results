% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_maf_triplet.R
\name{annotate_maf_triplet}
\alias{annotate_maf_triplet}
\title{Annotate MAF with triplet context}
\usage{
annotate_maf_triplet(
  maf,
  all_SNVs = TRUE,
  ref,
  alt,
  projection = "grch37",
  fastaPath,
  bsgenome_name,
  pyrimidine_collapse = FALSE
)
}
\arguments{
\item{maf}{MAF file (required columns: Reference_Allele, Tumor_Seq_Allele2)}

\item{all_SNVs}{To give us all the triplet sequences of SNVs and not
specifying any specific ref and alt alleles (default is TRUE)}

\item{ref}{Reference allele}

\item{alt}{Alternative allele}

\item{projection}{The genome build projection for the variants you are
working with (default is grch37)}

\item{fastaPath}{Can be a path to a FASTA file on a disk. When on GSC,
this is first attempted to be inferred from the gambl reference through
path specified in config. Local files are also accepted as value here.}

\item{bsgenome_name}{Name of a BSgenome data package (It has 4 or 5 parts,
separated by dot: 1st(BSgenome) . 2nd:name of organism(Hsapiens) .
3rd:name of genome provider (UCSC, NCBI, TAIR,...) . 4th:name of NCBI
assembly (e.g. GRCh38) or UCSC genome (e.g. hg38) . 5th(optional): If
the package contains masked sequences (masked))}

\item{pyrimidine_collapse}{Estimate mutation_strand and}
}
\value{
A data frame with two to three extra columns, in case
pyrimidine_collapse = FALSE, it will add triple sequence (seq) and the
strand (mutation_strand). When pyrimidine_collapse = T, another column
also will be added to the those extra columns which shows the mutation
(mutation)
}
\description{
Give triple sequence of mutated base with its adjacent bases (-1
and +1)
}
\details{
It gives the reference and alternative alleles and looks for the
rows of the data frame based on these values for + strand genes and their
complement alleles rows for - strand genes, then it can look for the adjacent
bases in that mutation position. Also, it can look for all the SNVs in
the MAF data frame and provide triple sequences for them (reverse complement
sequence for the - strand).
}
\examples{
maf <- GAMBLR.data::get_coding_ssm(these_sample_id = "DOHH-2")

annotate_maf_triplet(maf)
annotate_maf_triplet(maf, all_SNVs = FALSE, "C", "T")
annotate_maf_triplet(maf, ref = "C", alt = "T", pyrimidine_collapse = TRUE)

}
