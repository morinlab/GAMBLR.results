#! /usr/bin/env Rscript
# Usage: Rscript /path/to/GAMBLR.results/tools/logExampleOutputs.R
# If using renv, launch from your project root directory where renv.lock is located
# Do NOT run interactively

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_path <- normalizePath(script_path)

if(length(script_path) != 1){
	stop("script_path cannot be determined. Use Rscript /path/to/GAMBLR.results/tools/logExampleOutputs.R. Do not run interactively.")
}

setwd(dirname(dirname(script_path)))

log_file = "GAMBLR_examples_output.log"
options(width=2000)

sink(log_file)
print(paste("=== STARTED AT",Sys.time(),"==="))
devtools::run_examples()
print(paste("=== COMPLETED AT",Sys.time(),"==="))
sink()
