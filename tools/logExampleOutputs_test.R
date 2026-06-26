log_file = "GAMBLR_examples_output.log"
options(width = 2000)

sink(log_file, split = TRUE)
print(paste("=== STARTED AT", Sys.time(), "==="))

print("=== RUNNING devtools::run_examples() ===")
devtools::run_examples()

print(paste("=== COMPLETED AT", Sys.time(), "==="))
sink()

