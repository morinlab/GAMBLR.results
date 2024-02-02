## GABMLR.results - an R package with convenience functions for working with GAMBL results (restricted to useers with access to GSC).

## Installation

GAMBLR.results is an open-source package. It can be easily installed directly from GitHub:

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

## Running GAMBLR On Your Own Computer

If you don't have access to gphost on GSC, no worries, you can still execute GAMBLR functions in another way. Remote support was developed for this purpose. This section explains how to run GAMBLR remote on a *local machine* (i.e on your own computer). There are two different approaches to get this to work, both with its own advantages and limitations. We will be going over both in this next section.

### Approach 1 - Quick Start

This section details how to deploy GAMBLR with limited functionality. This approach requires either a working GSC VPN connection (or is directly accessible if connected to the GSC network).

#### Setup VPN Connection

1. You need a working GSC VPN connection to use this approach. For setting up a VPN connection see [this](https://www.bcgsc.ca/wiki/pages/viewpage.action?spaceKey=SysHelp&title=Learn+how+to+use+VPN) guide. Keep in mind that a **VPN connection is not needed** if your already connected to the GSC network.

#### Clone Repos, Update Paths, Install and Load R Packages

2. Clone [GAMBL](https://github.com/morinlab/gambl) and [GAMBLR](https://github.com/morinlab/GAMBLR) to your local computer. From your terminal run the following commands (folder structures can be whatever you want...)

```
mkdir ~/git_repos
cd ~/git_repos #set as working directiory
git clone https://github.com/morinlab/gambl
git clone https://github.com/morinlab/GAMBLR
```

3. Update the **paths** in your local [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) (GAMBLR) to point to the recently cloned, local **gambl** folder (repo_base). In your favorite text editor, edit the line shown below (under  *remote* ). Similarly, you will also need to edit the line above it to point to where you will eventually sync the GAMBL results.

```
remote:
    project_base: "/path/to/your/local/gambl_results_directory/"
    repo_base: "/path/to/your/local/gambl_repo/"
```

4. Set the **working directory** in Rstudio. Open Rstudio on your local machine and locate the repo you cloned previously.

```
setwd("~/git_repos/GAMBLR")
```

5. Install GAMBLR in your local R studio.

```
devtools::install()
```

6. Load packages.

```
library(GAMBLR)
```

#### Set Config To Remote

7. Execute the following in Rstudio console to make use of the updated paths in the [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) from step 3.

```
Sys.setenv(R_CONFIG_ACTIVE = "remote")
```

Alternatively, you can add the content of [~/git_repos/GAMBLR/.Rprofile](https://github.com/morinlab/GAMBLR/blob/master/.Rprofile) to your `~/.Rprofile` file. In this way, you do not need to enter the command above every time you start your R session (recommended).

```
cat ~/git_repos/GAMBLR/.Rprofile >> ~/.Rprofile
```

#### Run GAMBLR

8. Test if setup was successful (e.g call `get_gambl_metadata()` to retrieve meta data for all gambl samples).

```
get_gambl_metadata() %>%
  head()
```

### Approach 2 - The Full Installation (Snakemake)

This section details how to obtain GAMBLR with **full** functionality, using a dedicated snake file to retrieve all necessary files and dependencies.

#### Before You Get Started

1. Make sure you have a working SSH key  **setup with a pass phrase** . If not, follow instructions at [GSC Wiki](https://www.bcgsc.ca/wiki/login.action?os_destination=%2Fpages%2Fviewpage.action%3FspaceKey%3DSysHelp%26title%3DSetup%2BPassword-less%2BSSH&permissionViolation=true). Warning, this will **not** work with a pass phrase-less SSH connection.

#### Clone Repos and Set Up Environment

2. Clone [GAMBL](https://github.com/morinlab/gambl) and [GAMBLR](https://github.com/morinlab/GAMBLR).

```
mkdir ~/git_repos
cd ~/git_repos
git clone https://github.com/morinlab/gambl
git clone https://github.com/morinlab/GAMBLR
```

3. On your local machine, make a new directory called  **gambl_results** , for example.

```
mkdir ~/gambl_results/
```

4. In the [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) file of your local **GAMBLR** folder, update paths under the `remote` field to point to the recently cloned local **gambl** folder ( *repo_base* ) and recently created **gambl_results** folder ( *project_base* ). Also, update the *host* field to contain your username (you can use any gphost here). For example:

```
remote:
    project_base: "~/gambl_results/"
    repo_base: "~/git_repos/gambl/"
    ...
    host: "your_username@gphost01.bcgsc.ca"
```

5. Copy the following files (from your recently cloned [GAMBLR](https://github.com/morinlab/GAMBLR) directory) into the folder from the previous step; `config.yml` and `get_gambl_results.smk`.

```
cp ~/git_repos/GAMBLR/config.yml ~/gambl_results/
cp ~/git_repos/GAMBLR/get_gambl_results.smk ~/gambl_results/
```

6. Add ENVVARS bash/zsh environment variables to your bashrc/zsh or some other way that will ensure they're in your session (e.g. you can set them manually each time if you want, just make sure they are set). For example in your local terminal run the following commands (with updated values...).

```
export GSC_USERNAME="your_gsc_username"
export GSC_KEY="path_to_SSH_key_with_passphrase_from_step_1"
export GSC_PASSPHRASE="passpharase_from_step_1"
```

#### Install GAMBLR In Local Rstudio

7. Open **Rstudio** (locally) and set the working directory to the folder you downloaded in step 2 (in the Rstudio console) and install GAMBLR.

```
setwd("~/git_repos/GAMBLR")
```

8. Install and load GAMBLR into your local R session.

```
devtools::install()
```

#### Create and Setup Snakemake Environment

9. In the terminal on your local machine, create a new snakemake environment from the [get_gambl_results.yml](https://github.com/morinlab/GAMBLR/blob/master/get_gambl_results.yml) file ([get_gambl_results_linux.yml](https://github.com/morinlab/GAMBLR/blob/master/get_gambl_results_linux.yml) for Linux). Note that you can name this new environment whatever you would like. In this example, the new environment is called  **snakemake_gambl** .

```
cd ~/gambl_results
conda env create --name snakemake_gambl --file ~/git_repos/GAMBLR/get_gambl_results.yml
```

10. Activate this newly created snakemake environment with:

```
conda activate snakemake_gambl
```

11. Retrieve necessary files (download a local copy of all files needed to run a collection of GAMBLR functions). It's strongly advised to use `--cores 1` for this, since it seems to be the more stable option. In addition, if your sync gets interrupted, you only need restart the syncing of 1 file, compared to if you run on multiple cores.

```
snakemake -s get_gambl_results.smk --cores 1
```

#### Use GAMBLR Functions Locally

12. In Rstudio (local), open [test_remote.R](https://github.com/morinlab/GAMBLR/blob/master/resources/test_remote.R) in GAMBLR master folder.
13. Execute the following in Rstudio console to make use of the updated paths in the [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) from step 5.

```
Sys.setenv(R_CONFIG_ACTIVE = "remote")
```

Alternatively, you can add the content of [~/git_repos/GAMBLR/.Rprofile](https://github.com/morinlab/GAMBLR/blob/master/.Rprofile) to your `~/.Rprofile` file. In this way, you do not need to enter the command above every time you start your R session (recommended).

```
cat ~/git_repos/GAMBLR/.Rprofile >> ~/.Rprofile
```

14. Check what files (if any) are currently missing.

```
check_gamblr_config()
```

15. You should now be all set to explore a collection of GAMBLR function remotely on your local machine. For example you could try the following test code to ensure your setup was successful. For a set of comprehensive examples and tutorials, please refer to the [test_remote.R](https://github.com/morinlab/GAMBLR/blob/master/resources/test_remote.R) script.

```
get_gambl_metadata() %>%
  head()
```

 **Note** , if your seeing the following message when trying to use GAMBLR, please ensure that the config/gambl repo is set up properly (step  **5** ) and that you are using the *remote* mode (step  **13** ).

```
get_gambl_metadata(seq_type_filter = "capture") %>%
  pull(cohort) %>%
  table()

Error: '/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/metadata/gambl_all_outcomes.tsv' does not exist.
```
