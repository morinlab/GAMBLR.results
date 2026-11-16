## GABMLR.results - an R package for obtaining and manipulating GAMBL results (restricted to Morin lab members and collaborators).

## Installation

GAMBLR.results is an open-source package. It can be easily installed directly from GitHub using the command shown below. Please note that the functionality provided by this package will only work for users within the BC Cancer network who have been approved to access the underlying data sets. If this does not describe you, we recommend GAMBLR.open.

```
devtools::install_github("morinlab/GAMBLR.results", repos = BiocManager::repositories())
```

This will install the [GAMBLR.results](https://github.com/morinlab/GAMBLR.results) package with all necessary dependencies. It requires access to the GSC resources and is not intended to be used outside of GSC. If you are interested in standalone functionality, please refer to the documentation of the [GAMBLR.data](https://github.com/morinlab/GAMBLR.data) package or any other individual child package.

## Reproducing results from before the 2026-11 blacklist audit

In November 2026, positions shown to be real hotspots (e.g. FOXC1 S82N, HLA-B Y108*, CD70 splice sites) were taken off the SLMS-3 clean blacklists, and the affected MAFs were regenerated. The current results therefore include calls that earlier versions dropped. For results consistent with analyses done before the audit, pass `blacklist_version = "2022"` to any of the SSM getters (`get_coding_ssm`, `get_all_coding_ssm`, `get_ssm_by_sample(s)`, `get_ssm_by_region(s)`, `get_ssm_by_genes`), or apply it to calls you already have:

```
meta = get_gambl_metadata() %>% dplyr::filter(seq_type == "capture")
pre_audit = get_coding_ssm(these_samples_metadata = meta, this_seq_type = "capture",
                           blacklist_version = "2022")
# or, on an existing data frame of calls
pre_audit = apply_blacklist_version(current_calls, meta, blacklist_version = "2022")
```

The default, `"2026-11"`, returns the current calls. A call is dropped only if its sample's native-build blacklists (metadata `genome_build`) dropped that position, just as the pipeline did. See `?apply_blacklist_version` for details. The removed positions are read from `level_3/{grch37,hg38}_removed_from_clean_blacklist.tsv` (config `resources$blacklist$removed$template`), so sync those files if you work from a local copy of the results.

## Contributing

If you have access to gphost, the easiest way to obtain and contribute to GAMBLR is to do this via cloning the repository

```
cd
git clone git@github.com:morinlab/GAMBLR.results.git
```

In your R editor of choice, set your working directory to the place you just cloned the repo.

```
setwd("~/GAMBLR.results")
```

Install the package in R by running the following command (requires the devtools package)

```
devtools::install()
```

As GAMBL users (GAMBLRs, so to speak) rely on the functionality of this package, the Master branch is protected. All commits must be submitted via pull request on a branch. Please refer to the [GAMBL](https://github.com/morinlab/gambl#contribution-guidelines) documentation for details on how to do this.

