% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collate_results.R
\name{collate_results}
\alias{collate_results}
\title{Collate Results}
\usage{
collate_results(
  these_samples_metadata,
  write_to_file = FALSE,
  join_with_full_metadata = FALSE,
  sample_table,
  case_set,
  sbs_manipulation = "",
  seq_type_filter = "genome",
  from_cache = TRUE
)
}
\arguments{
\item{these_samples_metadata}{Optional argument to use a user specified metadata df, must include seq_type column to specify desired seq types. If not provided, defaults to all genome samples with \code{get_gambl_metadata}.}

\item{write_to_file}{Boolean statement that outputs tsv file (/projects/nhl_meta_analysis_scratch/gambl/results_local/shared/gambl_{seq_type_filter}_results.tsv) if TRUE, default is FALSE.}

\item{join_with_full_metadata}{Join with all columns of metadata, default is FALSE.}

\item{sample_table}{A data frame with sample_id as the first column (deprecated, use these_samples_metadata instead).}

\item{case_set}{Optional short name for a pre-defined set of cases.}

\item{sbs_manipulation}{Optional variable for transforming sbs values (e.g log, scale).}

\item{seq_type_filter}{Filtering criteria, default is genomes (deprecated. Include a \code{seq_type} column in these_samples_metadata to specify seq type).}

\item{from_cache}{Boolean variable for using cached results (/projects/nhl_meta_analysis_scratch/gambl/results_local/shared/gambl_{seq_type}_results.tsv), default is TRUE. If write_to_file is TRUE, this parameter auto-defaults to FALSE.}
}
\value{
A table keyed on biopsy_id that contains a bunch of per-sample results from GAMBL
}
\description{
Bring together all derived sample-level results from many GAMBL pipelines.
}
\details{
This function takes a data frame with sample IDs and seq types (in sample_id and seq_type columns) with the \code{these_samples_metadata} parameter and adds sample-level results from many of the available GAMBL pipelines.
Optional parameters are \code{join_with_full_metadata}. If \code{join_with_full_metadata} is set to TRUE, the function will return all metadata with \code{get_gambl_metadata}, allowing the user to extend the available information in a metadata table.
This function has also been designed so that it can get cached results, meaning that not all individual collate helper functions would have to be run to get results back.
To do so, run this function with \code{from_cache = TRUE} (default). In addition, it's also possible to regenerate the cached results, this is done by setting \code{write_to_file = TRUE},
This operation auto defaults \code{from_cache = FALSE}. \code{case_set} is an optional parameter available for subsetting the return to an already defined set of cases.
If a dataframe is not provided, the function will default to all genome metadata returned with \code{get_gambl_metadata}. For more information on how to get the most out of this function,
refer to function examples, vignettes and parameter descriptions.
}
\examples{
\dontrun{
#get collated results for all capture samples, using cached results
capture_metadata = get_gambl_metadata(dna_seq_type_priority = "capture") \%>\% dplyr::filter(seq_type == "capture")
capture_collated_everything = collate_results(these_samples_metadata = capture_metadata,
                                              from_cache = TRUE,
                                              write_to_file = FALSE)

#use an already subset metadata table for getting collated results (cached)
my_metadata = get_gambl_metadata()
fl_metadata = dplyr::filter(my_metadata, pathology == "FL")

fl_collated = collate_results(these_samples_metadata = fl_metadata,
                              write_to_file = FALSE,
                              from_cache = TRUE)
}
#use an already subset metadata table for getting collated results (without using cached results)
my_metadata = get_gambl_metadata()
fl_metadata = dplyr::filter(my_metadata, pathology == "FL")

fl_collated = collate_results(these_samples_metadata = fl_metadata,
                              write_to_file = FALSE,
                              from_cache = FALSE)
dplyr::select(fl_collated, 1:14) \%>\% head()
}
