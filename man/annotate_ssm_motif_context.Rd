% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_ssm_motif_context.R
\name{annotate_ssm_motif_context}
\alias{annotate_ssm_motif_context}
\title{Annotate mutations target motif}
\usage{
annotate_ssm_motif_context(
  maf,
  motif = "WRCY",
  index = 3,
  projection = "grch37",
  fastaPath
)
}
\arguments{
\item{maf}{MAF data frame (required columns: Reference_Allele, Chromosome, Start_Position, End_Position)}

\item{motif}{The motif sequence (default is WRCY)}

\item{index}{Position of the mutated allele in the motif}

\item{projection}{The genome build projection for the variants you are working with (default is grch37)}

\item{fastaPath}{Can be a path to a FASTA file}
}
\value{
A data frame with two extra columns (seq and motif).
}
\description{
Checks for the presence of mutations at a given motif
}
\details{
In positions that reference allele has been mutated, it will capture (motif length - 1) before and (motif length + 1) alleles after the mutated position.
Then, it looks for the presence of motif in the captured sequence and check if the mutation has occurred in the indexed position, it will return
SITE and if the the motif is present, but the mutation is not in the indexed position, it will return MOTIF. In other cases, it will return FALSE.
NA will be shown if the mutation is an indel mutation.
}
\examples{
my_maf <- get_coding_ssm()
annotate_ssm_motif_context(maf = my_maf, motif = "WRCY", index = 3)

}
