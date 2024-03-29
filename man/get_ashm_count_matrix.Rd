% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ashm_count_matrix.R
\name{get_ashm_count_matrix}
\alias{get_ashm_count_matrix}
\title{Get ASHM Count Matrix.}
\usage{
get_ashm_count_matrix(
  regions_bed,
  maf_data,
  these_samples_metadata,
  this_seq_type
)
}
\arguments{
\item{regions_bed}{A bed file with one row for each region.}

\item{maf_data}{Optionally provide a data frame in the MAF format, otherwise
the database will be used.}

\item{these_samples_metadata}{This is used to complete your matrix. All GAMBL
samples will be used by default. Provide a data frame with at least
sample_id for all samples if you are using non-GAMBL data.}

\item{this_seq_type}{The seq type to return results for. Only used if no
metadata is provided with these_samples_metadata.}
}
\value{
matrix
}
\description{
Prepare a matrix with one row per sample and one column per
region using a set of hypermutated regions.
}
\details{
Values are the number of mutations in that patient in the region.
}
\examples{
regions_bed <- dplyr::mutate(
     GAMBLR.data::grch37_ashm_regions,
     name = paste(gene, region, sep = "_")
)

matrix <- get_ashm_count_matrix(
     regions_bed = regions_bed,
     this_seq_type = "genome"
)

}
