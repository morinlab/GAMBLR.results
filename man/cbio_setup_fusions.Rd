% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cbio_setup_fusions.R
\name{cbio_setup_fusions}
\alias{cbio_setup_fusions}
\title{Setup Fusions (cBioPortal).}
\usage{
cbio_setup_fusions(
  short_name = "GAMBL",
  human_friendly_name = "GAMBL data",
  project_name = "gambl_genome",
  description = "GAMBL data from genome",
  gambl_maf = "maf_slms3_hg19",
  gambl_icgc_maf = "maf_slms3_hg19_icgc",
  out_dir
)
}
\arguments{
\item{short_name}{A concise name for your portal project.}

\item{human_friendly_name}{A slightly more verbose name for your project.}

\item{project_name}{Unique ID for your project.}

\item{description}{A verbose description of your data set.}

\item{gambl_maf}{maf origin.}

\item{gambl_icgc_maf}{ICGC maf origin.}

\item{out_dir}{The full path to the base directory where the files are being
created.}
}
\value{
A vector of sample_id for the patients that have been included.
}
\description{
Annotate SVs and create the input for fusions to be displayed in
cBioPortal instance.
}
\details{
This function calls \link{get_combined_sv} and runs
annotation on the returned data frame. Should be run as the next step after
running \link{cbio_setup_study}. Note that the parameters called
with this function has to match the parameter arguments of
\link{cbio_setup_study}, i.e if \code{short_name} is for
\link{cbio_setup_study} is "GAMBL", then the \code{short_name} in
\link{cbio_setup_fusions} also has to be "GAMBL", etc.
}
\examples{
\dontrun{
fusion_ids = cbio_setup_fusions(out_dir = "GAMBLR/cBioPortal/instance01/")
}
}
\keyword{internal}
