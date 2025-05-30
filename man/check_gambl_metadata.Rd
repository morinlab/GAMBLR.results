% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_gambl_metadata.R
\name{check_gambl_metadata}
\alias{check_gambl_metadata}
\title{Check GAMBL Metadata.}
\usage{
check_gambl_metadata(
  metadata_df,
  to_check = "all",
  show_details = FALSE,
  fix = FALSE
)
}
\arguments{
\item{metadata_df}{Data frame output by \code{get_gambl_metadata} or some other source of metadata you plan to use.}

\item{to_check}{Specify one of "uniqueness", "colours" or "completeness" or leave empty to check all.}

\item{show_details}{Set to TRUE if you want the gory details about issues that are identified.}

\item{fix}{After identifying an issue, rerun this function with fix=TRUE to address errors (when possible). Currently this doesn't do anything. That's how I roll.}
}
\description{
Check GAMBL or other metadata for compatibility with various features.
}
\details{
Give this function a metadata output, preferably from \link{get_gambl_metadata}, or any other source.
The function then checks for duplicate sample IDs, colours for all values in all columns that map to a colour with \code{map_metadata_to_colours} and missing columns.
}
\examples{
\dontrun{
this_metadata = get_gambl_metadata()
check_gambl_metadata(metadata_df = this_metadata)
}
}
\keyword{internal}
