# Prefer the source tree's inst/extdata so the tests run in CI without
# installing the package; fall back to the installed copy under R CMD check
extdata_file <- function(name) {
  src <- test_path("..", "..", "inst", "extdata", name)
  if (file.exists(src)) src else system.file("extdata", name, package = "GAMBLR.results")
}

# Keys the package reads, as listed in the manifest GAMBLR checks on attach
required_keys <- function() {
  readLines(extdata_file("required_config_keys.txt"))
}

has_config_key <- function(config_list, key) {
  x <- config_list
  for (part in strsplit(key, "$", fixed = TRUE)[[1]]) {
    if (!is.list(x) || !part %in% names(x)) {
      return(FALSE)
    }
    x <- x[[part]]
  }
  TRUE
}

test_that("template config.yml has every key in required_config_keys.txt", {
  template <- extdata_file("config.yml")
  expect_true(nzchar(template))

  config_list <- config::get(file = template, config = "default")
  keys <- required_keys()
  missing_keys <- keys[!vapply(keys, has_config_key, logical(1), config_list = config_list)]
  expect_equal(missing_keys, character(0))
})

test_that("required_config_keys.txt matches the keys read in R/", {
  # Needs the package source, so only runs under devtools::test(), not R CMD check
  pkg_root <- test_path("..", "..")
  skip_if_not(dir.exists(file.path(pkg_root, "R")), "package source not available")
  skip_if_not(file.exists(file.path(pkg_root, "tools", "update_required_config_keys.R")))

  tools_env <- new.env()
  sys.source(file.path(pkg_root, "tools", "update_required_config_keys.R"), envir = tools_env)
  current <- tools_env$find_required_config_keys(file.path(pkg_root, "R"))
  manifest <- readLines(file.path(pkg_root, "inst", "extdata", "required_config_keys.txt"))

  not_in_manifest <- setdiff(current, manifest)
  no_longer_read <- setdiff(manifest, current)
  expect(
    length(not_in_manifest) == 0 && length(no_longer_read) == 0,
    paste0(
      "required_config_keys.txt is out of date. Run `Rscript tools/update_required_config_keys.R`.\n",
      "Not in manifest: ", paste(not_in_manifest, collapse = ", "), "\n",
      "No longer read: ", paste(no_longer_read, collapse = ", ")
    )
  )
})
