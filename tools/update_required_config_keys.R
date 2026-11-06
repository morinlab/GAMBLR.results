# Regenerate inst/extdata/required_config_keys.txt by scanning R/ for every
# literal config key this package reads. GAMBLR checks the user's config.yml
# against this manifest when it is attached, and tests/testthat/test-config.R
# fails if the manifest is out of date.
#
# Run from the GAMBLR.results package root after adding or renaming config keys:
#   Rscript tools/update_required_config_keys.R

# Keys that are only defined for a non-default config (e.g. remote) and
# should not be required of everyone
exclude_config_keys <- c("host")

find_required_config_keys <- function(r_dir = "R", exclude = exclude_config_keys) {
  src <- unlist(lapply(list.files(r_dir, pattern = "\\.R$", full.names = TRUE),
                       readLines, warn = FALSE))
  src <- gsub("[[:space:]]+", "", src)

  match_all <- function(pattern, text) {
    unlist(regmatches(text, gregexpr(pattern, text, perl = TRUE)))
  }

  # check_config_and_value("a$b$c")
  keys <- sub('.*\\("([^"]+)".*', "\\1",
              match_all('check_config_and_value\\("[^"]+"', src))

  # config::get("a")$b$c, skipping lookups pinned to a non-default config
  gets <- match_all('config::get\\("[A-Za-z_]+"(,config="[a-z]+")?\\)(\\$[A-Za-z0-9_]+)*', src)
  gets <- gets[!grepl('config="(?!default)', gets, perl = TRUE)]
  keys <- c(keys, paste0(sub('config::get\\("([A-Za-z_]+)".*', "\\1", gets),
                         sub('^[^)]*\\)', "", gets)))

  # get_template_wildcards("a") reads a, get_template_wildcards("a", x) reads a_wildcards
  wild <- match_all('get_template_wildcards\\("[A-Za-z_]+"[,)]', src)
  keys <- c(keys, paste0(sub('.*\\("([A-Za-z_]+)".*', "\\1", wild),
                         ifelse(endsWith(wild, ","), "_wildcards", "")))

  keys <- setdiff(unique(keys), exclude)

  # A leaf key implies its parents, so keep only the most specific paths
  is_parent <- vapply(keys, function(k) any(startsWith(keys, paste0(k, "$"))), logical(1))
  # radix sorts the same in every locale, so the manifest doesn't churn
  sort(keys[!is_parent], method = "radix")
}

# Write the manifest only when run as a script, not when sourced by the tests
if (sys.nframe() == 0) {
  keys <- find_required_config_keys()
  dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
  writeLines(keys, "inst/extdata/required_config_keys.txt")
  message("Wrote ", length(keys), " keys to inst/extdata/required_config_keys.txt")
}
