% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_all_coding_ssm.R
\name{get_all_coding_ssm}
\alias{get_all_coding_ssm}
\title{Get all Coding SSMs}
\usage{
get_all_coding_ssm(these_samples_metadata = NULL, include_silent = FALSE, ...)
}
\arguments{
\item{these_samples_metadata}{Supply a metadata table containing the sample/seq_type combinations you want.}

\item{include_silent}{If set to TRUE, silent/synonymous mutations in the coding regions will also be returned.}
}
\value{
A data frame containing all the MAF data columns (one row per mutation).
}
\description{
Retrieve all coding SSMs from GAMBL in MAF-like format, regardless of seq_type.
}
\details{
Effectively retrieve coding SSM calls from one or all DNA seq_type. For additional optional arguments, see \link{get_coding_ssm}
}
\examples{

maf_data_all_exome_and_genome = get_coding_ssm(these_samples_metadata = get_gambl_metadata(seq_type_filter=c("genome","capture")))


}
