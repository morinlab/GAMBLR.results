## GABMLR.results - access and manipulate GAMBL results in R

## Installation

GAMBLR.results is an open-source package. It can be easily installed directly from GitHub using the command shown below. Please note that the functionality provided by this package will only work for users within the BC Cancer network who have been approved to access the underlying data sets. If this does not describe you, we recommend GAMBLR.open.

```
devtools::install_github("morinlab/GAMBLR.results", repos = BiocManager::repositories())
```

This will install the [GAMBLR.results](https://github.com/morinlab/GAMBLR.results) package with all necessary dependencies. It requires access to the GSC resources and is not intended to be used outside of GSC. If you are interested in standalone functionality, please refer to the documentation of the [GAMBLR.data](https://github.com/morinlab/GAMBLR.data) package or any other individual child package.

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

