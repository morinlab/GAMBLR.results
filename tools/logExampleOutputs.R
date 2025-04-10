log_file = "GAMBLR_examples_output.log"
options(width=2000)

sink(log_file)
print(paste("=== STARTED AT",Sys.time(),"==="))
devtools::run_examples()
print(paste("=== COMPLETED AT",Sys.time(),"==="))
sink()
