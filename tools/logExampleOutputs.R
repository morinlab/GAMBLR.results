library(optparse)

option_list = list(
	make_option(c("-e", "--renv"), type="logical", default=FALSE, action="store_true", 
		help="Flag to specify whether to use renv or not. [default %default]", metavar="character"),
	make_option(c("-p", "--path"), type="character", default=NULL, action="store_true", 
		help="Full path to the directory with the renv files.", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

current_path <- getwd()

use_renv <- opt$renv
if(use_renv){
	if(is.null(opt$path)){
		stop("Pleaser provide the full path to the directory with the renv files.")
	}else if(!dir.exists(opt$path)){
		stop(paste("Directory,", opt$path, "does not exist."))
	}else{
		cat("Using renv\n")
		Sys.setenv(RENV_PROJECT = opt$path)
		setwd(opt$path)
		renv::load()
		setwd(current_path)
	}
}

log_file = "GAMBLR_examples_output.log"
options(width=2000)

sink(log_file)
print(paste("=== STARTED AT",Sys.time(),"==="))
devtools::run_examples()
print(paste("=== COMPLETED AT",Sys.time(),"==="))
sink()
