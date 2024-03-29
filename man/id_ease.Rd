% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/id_ease.R
\name{id_ease}
\alias{id_ease}
\alias{id_ease,}
\alias{id}
\alias{ease}
\title{ID Ease}
\usage{
id_ease(
  these_samples_metadata = NULL,
  these_sample_ids = NULL,
  this_seq_type = c("genome", "capture"),
  verbose = FALSE
)
}
\arguments{
\item{these_samples_metadata}{An optional data frame with metadata, subset to sample IDs of interest.
If not provided will retrieve GAMBL metadata for all available samples.}

\item{these_sample_ids}{Optional character vector of GAMBL sample IDs.}

\item{this_seq_type}{The seq type of interest. Default is both genome and exome, with priority for genome when a sample has >1 seq_type.}

\item{verbose}{Set to FALSE to limit the information that gets printed to the console. Default is FALSE.}
}
\value{
Metadata (data frame).
}
\description{
Internal convenience function that standardize the way GAMBLR functions deals with sample IDs (these_sample_ids)
and metadata (these_samples_metadata).
}
\details{
This function can take sample IDs as a vector of characters, or a metadata table in data frame format.
If no sample IDs are provided to the function, the function will operate on all gambl sample IDs available for the given seq type.
It is highly recommended to run this function with \code{verbose = TRUE}.
Since this will not only improve the overall logic on how the function operates.
But also might help with debugging functions that are internally calling this function.
The function also performs sanity checks and notifies the user if any of the requested sample IDs are not found in the metadata.
In addition, the function also notifies the dimensions of the returned object, providing further insight to what is returned.
As with all GAMBLR functions, providing a curated metadata table to any GAMBLR function (as opposed to a vector of IDs) is the safest way to ensure you get the expected result.
}
\examples{
#load packages
library(dplyr)

#give the function nothing (i.e return all sample IDs in the metadata for the default seq type)
#return metadata for all samples in the default seq type
all_meta = id_ease()

#return metadata based on a sample ID
sample_meta = id_ease(these_sample_ids = "94-15772_tumorA")

#return sample IDs based on an already filtered metadata
this_metadata = get_gambl_metadata(seq_type_filter = "genome") \%>\% 
  head(5)

these_ids = id_ease(these_samples_metadata = this_metadata)

}
