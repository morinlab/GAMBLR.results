% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_gene_expression.R
\name{check_gene_expression}
\alias{check_gene_expression}
\title{Check Gene Expression.}
\usage{
check_gene_expression(verbose = F, show_linkages = F, ...)
}
\arguments{
\item{verbose}{Set to TRUE mainly for debugging}

\item{show_linkages}{Set to TRUE to link every row to an available capture and genome sample
using get_gambl_metadata to prioritize each to at most one per biopsy_id}

\item{...}{Optional parameters to pass along to get_gambl_metadata (only used if show_linkages = TRUE)}
}
\value{
A data frame with a row for each non-redundant RNA-seq result and the following columns:

\describe{
\item{mrna_sample_id}{The unique sample_id value that will match a single row from the GAMBL metadata where seq_type is mrna. }
\item{biopsy_id}{The unique identifier for the source of nucleic acids.}
\item{sample_id}{Identical to mrna_sample_id}
\item{capture_sample_id}{When this biopsy has capture/exome data in the GAMBL metadata, the value will be the sample_id for that data. NA otherwise.}
\item{genome_sample_id}{When this biopsy has genome data in the GAMBL metadata, the value will be the sample_id for that data. NA otherwise.}
\item{patient_id}{The anonymized unique identifier for this patient. For BC samples, this will be Res ID.}
\item{seq_type}{The assay type used to produce this data (will always be "mrna" in this case)}
\item{protocol}{Specifies the RNA-seq library construction protocol.}
\item{ffpe_or_frozen}{Specifies the way the source of nucleic acids was preserved. Either FFPE or frozen.}
}
}
\description{
This function determines which samples have expression data available in the merge and drop redundant data while consistently prioritizing by protocol and nucleic acid source.
}
