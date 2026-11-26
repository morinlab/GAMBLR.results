# A MAF path or its gzipped copy: `path` if that file exists, otherwise `path`.gz if that exists, otherwise `path`
# unchanged (so a caller's own missing-file handling still reports the expected name). Lets local copies of MAFs
# (e.g. rarely used per-sample MAFs) be gzipped without breaking anything that reads them: data.table::fread()
# and readr (and so GAMBLR.utils::fread_maf()) read .gz files transparently.
resolve_maf_path = function(path){
  if(length(path) != 1 || is.na(path) || file.exists(path)) return(path)
  gz = paste0(path, ".gz")
  if(file.exists(gz)) gz else path
}

# TRUE if the MAF at `path`, or its gzipped copy, exists
maf_file_exists = function(path){
  file.exists(resolve_maf_path(path))
}

# shell command that streams a MAF's lines through `filter_cmd` (e.g. a grep), decompressing it first if gzipped
maf_stream_cmd = function(path, filter_cmd){
  if(grepl("\\.gz$", path)){
    paste("gzip -dc", shQuote(path, type = "sh"), "|", filter_cmd)
  }else{
    paste(filter_cmd, shQuote(path, type = "sh"))
  }
}
