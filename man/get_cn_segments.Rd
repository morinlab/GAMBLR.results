% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cn_segments.R
\name{get_cn_segments}
\alias{get_cn_segments}
\title{Get CN Segments.}
\usage{
get_cn_segments(
  region,
  chromosome,
  qstart,
  qend,
  projection = "grch37",
  this_seq_type = "genome",
  with_chr_prefix = FALSE,
  streamlined = FALSE,
  from_flatfile = TRUE
)
}
\arguments{
\item{region}{Region formatted like chrX:1234-5678 or X:1234-56789.}

\item{chromosome}{The chromosome you are restricting to. Required parameter if region is not specified.}

\item{qstart}{Start coordinate of the range you are restricting to. Required parameter if region is not specified.}

\item{qend}{End coordinate of the range you are restricting to. Required parameter if region is not specified.}

\item{projection}{Selected genome projection for returned CN segments. Default is "grch37".}

\item{this_seq_type}{Seq type for returned CN segments. One of "genome" (default) or "capture".}

\item{with_chr_prefix}{Boolean parameter for toggling if chr prefixes should be present in the return, default is FALSE.}

\item{streamlined}{Return a basic rather than full MAF format. Default is FALSE.}

\item{from_flatfile}{Set to TRUE by default.}
}
\value{
A data frame with CN segments for the specified region.
}
\description{
Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range.
}
\details{
This function returns CN segments for s specified region.
There are multiple ways a region can be specified.
For example, the user can provide the full region in a "region" format (chr:start-end) to the \code{region} parameter.
Or, the user can provide chromosome, start and end coordinates individually with \code{chr}, \code{start}, and \code{end} parameters.
For more usage examples, refer to the parameter descriptions and examples in the vignettes.
Is this function not what you are looking for? Try one of the following, similar, functions; \link{assign_cn_to_ssm}, \link{get_cn_states}, \link{get_sample_cn_segments}
}
\examples{
#Example using chromosome, qstart and qend parameters:
segments_region_grch37 = get_cn_segments(chromosome = "chr8",
                                         qstart = 128723128,
                                         qend = 128774067)
# Example for the capture samples:
capture_segments_region_grch37 = get_cn_segments(
 chromosome = "chr8",
 qstart = 128723128,
 qend = 128774067,
 this_seq_type = "capture"
)

#Example using the regions parameter:
segments_region_hg38 = get_cn_segments(region = "chr8:128,723,128-128,774,067",
                                       projection = "hg38",
                                       with_chr_prefix = TRUE)

}
