% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_hotspots.R
\name{annotate_hotspots}
\alias{annotate_hotspots}
\title{Annotate Hotspots.}
\usage{
annotate_hotspots(
  mutation_maf,
  recurrence_min = 5,
  analysis_base = c("FL--DLBCL", "BL--DLBCL"),
  p_thresh = 0.05
)
}
\arguments{
\item{mutation_maf}{A data frame in MAF format.}

\item{recurrence_min}{minimum number of recurrences for mutation to be included, default is 5.}

\item{analysis_base}{Base name for hot spot output directory.}

\item{p_thresh}{P value threshold, default is 0.05.}
}
\value{
The same data frame with one additional column "hot_spot".
}
\description{
Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations.
}
\details{
This function takes an already loaded MAF data frame with the \code{mutation_maf} parameter.
The user can then control the minimum number of recurrences for mutations to be included with \code{recurrance_min},
The default is 5. \code{analysis_base} controls the base name go hotspot output directory.
Lastly, \code{p_thresh} sets the p value threshold, default is 0.05.
}
\examples{
my_metadata = suppressMessages(get_gambl_metadata())
# get a few SSMs to annotate
some_coding_ssm = get_coding_ssm(these_samples_metadata = my_metadata,
                                projection = "grch37",
                                this_seq_type = "genome") \%>\% 
                  dplyr::filter(Hugo_Symbol \%in\% c("EZH2","MEF2B","MYD88","KMT2D")) \%>\%
                  dplyr::arrange(Hugo_Symbol)
# peek at the data
dplyr::select(some_coding_ssm,1:10,37) \%>\% head()

hot_ssms = annotate_hotspots(some_coding_ssm)
hot_ssms \%>\% 
   dplyr::filter(!is.na(hot_spot)) \%>\% 
   dplyr::select(1:10,37,hot_spot) 

\dontrun{
#This example will raise an error due to the user supplying an unsupported genome build:
more_coding_ssm = get_coding_ssm(
                                these_samples_metadata = my_metadata,
                                projection = "hg38",
                                this_seq_type = "capture") \%>\% 
                  dplyr::filter(Hugo_Symbol \%in\% c("EZH2","MEF2B","MYD88","KMT2D")) \%>\%
                  dplyr::arrange(Hugo_Symbol)
# peek at the data
dplyr::select(more_coding_ssm,1:10,37) \%>\% head()

more_hot_ssms = annotate_hotspots(more_coding_ssm)
more_hot_ssms \%>\% 
   dplyr::filter(!is.na(hot_spot)) \%>\% 
   dplyr::select(1:10,37,hot_spot) 
}
}
