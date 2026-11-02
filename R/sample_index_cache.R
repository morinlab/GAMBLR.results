#' @title Get the samples present in a merged results file
#'
#' @description Return the sample IDs present in a merged results file (e.g.
#' the indexed SSM MAF, or the combined SV file), from a small index file
#' that is built once and rebuilt only when the results file changes. A
#' region query of such a file can't tell a sample with no events in the
#' region from a sample that isn't in the file at all; this index can.
#'
#' @details The index is a text file with one sample ID per line, named after
#' the results file (`<file>.samples.txt`): next to it when that directory is
#' writable, otherwise in the per-user cache directory
#' (`tools::R_user_dir("GAMBLR", "cache")`). It is (re)built when missing or
#' older than the results file. Building it reads the whole file once, in a
#' single streaming pass that keeps only the set of IDs in memory (for
#' compressed files this uses `gzip` and `awk`; several minutes for the
#' genome MAF). Lines starting with "#" and the header line are skipped.
#'
#' @param file Path to the results file (plain or gzip/bgzip compressed).
#' @param sample_column Name of the column with the sample IDs. Default
#' "Tumor_Sample_Barcode".
#' @param refresh Set to TRUE to rebuild the index.
#' @param verbose Report where the index is read from or written to.
#'
#' @return A character vector of sample IDs.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' maf = paste0(check_config_and_value("project_base"),
#'              glue::glue(check_config_and_value("results_flatfiles$ssm$template$merged$augmented"),
#'                         seq_type = "genome", projection = "grch37"), ".bgz")
#' length(get_file_samples(maf))
#' }
get_file_samples = function(file, sample_column = "Tumor_Sample_Barcode",
                            refresh = FALSE, verbose = FALSE){
  if(!file.exists(file)) stop("file not found: ", file)
  shared = paste0(file, ".samples.txt")
  user_dir = tools::R_user_dir("GAMBLR", "cache")
  user = file.path(user_dir, paste0(gsub("[^A-Za-z0-9._-]", "_", file), ".samples.txt"))
  fresh = function(f) file.exists(f) && file.mtime(f) >= file.mtime(file)
  if(!refresh){
    for(f in c(shared, user)){
      if(fresh(f)){
        if(verbose) message("sample index from cache: ", f)
        return(readLines(f))
      }
    }
  }
  message("building the sample index for ", basename(file),
          " (one-off; can take several minutes for large files)")
  ids = .scan_file_samples(file, sample_column)
  target = if(file.access(dirname(file), 2) == 0) shared else {
    dir.create(user_dir, recursive = TRUE, showWarnings = FALSE); user
  }
  ok = tryCatch({ writeLines(ids, target); TRUE }, error = function(e) FALSE)
  if(verbose) message(if(ok) paste("sample index cached in", target) else "sample index could not be cached")
  ids
}

# Unique values of one column of a (possibly compressed) tab-separated file,
# in one streaming pass that holds only the set of values in memory
.scan_file_samples = function(file, sample_column){
  compressed = grepl("\\.(gz|bgz)$", file)
  con = if(compressed) gzfile(file, "r") else file(file, "r")
  header = NULL
  repeat{
    line = readLines(con, n = 1)
    if(length(line) == 0) break
    if(!startsWith(line, "#")){ header = strsplit(line, "\t")[[1]]; break }
  }
  close(con)
  col = match(sample_column, header)
  if(is.na(col)) stop("no column ", sample_column, " in ", file)
  reader = if(compressed) paste("gzip -dc", shQuote(file)) else paste("cat", shQuote(file))
  awk = paste0("awk -F'\\t' 'substr($0,1,1) != \"#\" && $", col, " != \"", sample_column,
               "\" && !($", col, " in s) {s[$", col, "]=1; print $", col, "}'")
  ids = tryCatch(system(paste(reader, "|", awk), intern = TRUE),
                 error = function(e) NULL, warning = function(w) NULL)
  if(is.null(ids)){
    # no shell tools: read the column in R (slower, more memory)
    ids = unique(as.character(readr::read_tsv(file, col_select = dplyr::all_of(sample_column),
                                              comment = "#", lazy = TRUE, progress = FALSE)[[1]]))
  }
  sort(unique(ids[nzchar(ids)]))
}
