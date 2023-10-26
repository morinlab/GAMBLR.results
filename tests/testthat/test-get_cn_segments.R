#load packages
library(testthat)

test_that("Returned rows an columns consitency check", {
  expect_equal(nrow(get_cn_segments(region = "chr8:128,723,128-128,774,067")), 1745) #grch37
  expect_equal(nrow(get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38")), 1751) #hg38
  expect_equal(nrow(get_cn_segments(region = "chr8:128,723,128-128,774,067", this_seq_type = "capture")), 3024) #grch37
  expect_equal(nrow(get_cn_segments(region = "chr8:128,723,128-128,774,067", this_seq_type = "capture", projection = "hg38")), 3082) #hg38
  expect_equal(nrow(get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38")), 1751) #hg38
  expect_equal(ncol(get_cn_segments(region = "chr8:128,723,128-128,774,067")), 7) #grch37
  expect_equal(ncol(get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38")), 7) #hg38
})


test_that("Compare the return for region parameter vs chromosome, qstart and qend + different region formats", {
  expect_identical(get_cn_segments(chromosome = "chr8", qstart = 128723128, qend = 128774067),
                   get_cn_segments(region = "chr8:128,723,128-128,774,067"),
                   get_cn_segments(region = "chr8:128723128-128774067"),
                   get_cn_segments(region = "8:128723128-128774067"))
  
  expect_identical(get_cn_segments(chromosome = "chr8", qstart = 128723128, qend = 128774067, this_seq_type = "capture", projection = "hg38"),
                   get_cn_segments(region = "chr8:128,723,128-128,774,067", this_seq_type = "capture", projection = "hg38"),
                   get_cn_segments(region = "chr8:128723128-128774067", this_seq_type = "capture", projection = "hg38"),
                   get_cn_segments(region = "8:128723128-128774067", this_seq_type = "capture", projection = "hg38"))
})


test_that("Chr prefixes only show up where expected", {
  expect_true(all(grepl("chr", get_cn_segments(region = "chr8:128,723,128-128,774,067", with_chr_prefix = TRUE)[,"chrom"]))) #using grch37 (default), should be easy
  expect_true(all(grepl("chr", get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38", with_chr_prefix = TRUE)[,"chrom"]))) #what if we request projection that is already chr prefixed, will we get double chr prefixes?
  expect_true(all(!grepl("chr", get_cn_segments(region = "chr8:128,723,128-128,774,067", with_chr_prefix = FALSE)[,"chrom"]))) #request a projection that is not prefixed, and keep with_chr_prefixes = FALSE
  expect_true(all(!grepl("chr", get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "hg38", with_chr_prefix = FALSE)[,"chrom"]))) #remove chr prefixes from hg38 return
})


test_that("Is the streamlined option only returning the two expected columns", {
  expect_equal(ncol(get_cn_segments(region = "8:128723128-128774067", streamlined = TRUE)), 2)
  expect_true(all(c("ID", "CN") %in% names(get_cn_segments(region = "8:128723128-128774067", streamlined = TRUE))))
  expect_equal(ncol(get_cn_segments(region = "8:128723128-128774067", streamlined = TRUE, projection = "hg38")), 2)
  expect_true(all(c("ID", "CN") %in% names(get_cn_segments(region = "8:128723128-128774067", streamlined = TRUE, projection = "hg38"))))  
})
