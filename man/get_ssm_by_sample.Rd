% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ssm_by_sample.R
\name{get_ssm_by_sample}
\alias{get_ssm_by_sample}
\title{Get SSM By Sample.}
\usage{
get_ssm_by_sample(
  these_samples_metadata,
  tool_name = "slms_3",
  projection = "grch37",
  augmented = TRUE,
  flavour = "clustered",
  min_read_support = 3,
  basic_columns = TRUE,
  maf_cols = NULL,
  verbose = FALSE,
  this_sample_id,
  this_seq_type
)
}
\arguments{
\item{these_samples_metadata}{Required. A single row of metadata specifying
which sample_id and seq_type you desire the mutations from}

\item{tool_name}{The name of the variant calling pipeline (currently
only slms_3 is supported).}

\item{projection}{The projection genome build. Supports hg38 and grch37}

\item{augmented}{default: TRUE. Set to FALSE if you instead want the original
MAF from each sample for multi-sample patients instead of the augmented MAF.}

\item{flavour}{Currently this function only supports one flavour option but
this feature is meant for eventual compatibility with additional variant
calling parameters/versions.}

\item{min_read_support}{Only returns variants with at least this many
reads in t_alt_count (for cleaning up augmented MAFs).}

\item{basic_columns}{Return first 43 columns of MAF rather than full
details. Default is TRUE.}

\item{maf_cols}{if basic_columns is set to FALSE, the user can specify
which columns to be returned within the MAF. This parameter can either
be a vector of indexes (integer) or a vector of characters.}

\item{verbose}{Enable for debugging/noisier output.}

\item{this_sample_id}{Optional. Can be used to subset a multi-row these_samples_metadata
table to a single row.}

\item{this_seq_type}{Deprecated. Inferred from these_samples_metadata}
}
\value{
data frame in MAF format.
}
\description{
Get the SSMs (i.e. load MAF) for a single sample.
}
\details{
This was implemented to allow flexibility because there are some
samples that we may want to use a different set of variants than those in
the main GAMBL merge.
The current use case is to allow a force_unmatched output to be used
to replace the SSMs from the merge for samples with known contamination
in the normal.
This will also be useful to apply a blacklist to individual MAFs when coupled
with \link{annotate_ssm_blacklist}.
Is this function not what you are looking for? Try one of the related
functions:
\link{get_coding_ssm}, \link{get_coding_ssm_status},
\link{get_ssm_by_samples},
\link{get_ssm_by_region}, \link{get_ssm_by_regions}
}
\examples{

maf_samp = GAMBLR.results:::get_ssm_by_sample(
  get_gambl_metadata() \%>\% dplyr::filter(sample_id=="13-27975_tumorA"),
  augmented = FALSE
)
nrow(maf_samp)
maf_samp_aug = GAMBLR.results:::get_ssm_by_sample(
  get_gambl_metadata() \%>\% dplyr::filter(sample_id=="13-27975_tumorA"),
  augmented = TRUE
)
nrow(maf_samp_aug)


 some_maf = GAMBLR.results:::get_ssm_by_sample(
                          these_samples_metadata = get_gambl_metadata() \%>\%
                            dplyr::filter(sample_id == "HTMCP-01-06-00485-01A-01D",
                                     seq_type == "genome"),
                         projection = "hg38")
 dplyr::select(some_maf,1:10)
 

}
