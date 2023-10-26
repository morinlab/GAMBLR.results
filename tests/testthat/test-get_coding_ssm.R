#load packages
library(testthat)

test_that("Check for rows and column consistencies", {
  expect_equal(nrow(suppressWarnings(get_coding_ssm())), 211893)
  expect_equal(nrow(suppressWarnings(get_coding_ssm(this_seq_type = "capture"))), 742418)
  expect_equal(nrow(suppressWarnings(get_coding_ssm(these_samples_metadata = get_gambl_metadata()))), 211893)
  expect_equal(ncol(suppressWarnings(get_coding_ssm())), 45)
  expect_equal(ncol(suppressWarnings(get_coding_ssm(basic_columns = FALSE))), 116)
  expect_equal(nrow(suppressWarnings(get_coding_ssm(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2")))), 234)
  expect_equal(nrow(suppressWarnings(get_coding_ssm(projection = "hg38"))), 213279)
  expect_equal(nrow(suppressWarnings(get_coding_ssm(projection = "hg38", these_samples_metadata = get_gambl_metadata()))), 213279)
  expect_equal(ncol(suppressWarnings(get_coding_ssm(projection = "hg38"))), 45)
  expect_equal(nrow(suppressWarnings(get_coding_ssm(projection = "hg38", these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2")))), 234)
})


test_that("Does inlcude_silent work as intended", {
  expect_true(any(grepl("Silent", get_coding_ssm(include_silent = TRUE)[,"Variant_Classification"]))) #include `Variant_Classification == Silent`
  expect_true(all(!grepl("Silent", get_coding_ssm(include_silent = FALSE)[,"Variant_Classification"]))) #exclude `Variant_Classification == Silent`
})


test_that("Can we specify what MAF columns we want back?", {
  expect_true(all(c("Variant_Type", "Chromosome") %in% names(suppressWarnings(get_coding_ssm(basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome"))))))
  expect_equal(ncol(suppressWarnings(get_coding_ssm(basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome")))), 2)
  expect_true(all(c("Variant_Type", "Chromosome") %in% names(suppressWarnings(get_coding_ssm(projection = "hg38", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome"))))))
  expect_equal(ncol(suppressWarnings(get_coding_ssm(projection = "hg38", basic_columns = FALSE, maf_cols = c("Variant_Type", "Chromosome")))), 2)
})


test_that("Do the min vaf filters work as advertised?", {
  expect_gte(min(get_coding_ssm(min_read_support = 100)[,"t_alt_count"]), 100)
})


test_that("Do the built in metadata subset options work as advertised?", {
  #exclude DLBCL_cell_line samples (cohort)
  dlbcl_cell_lines_samples = get_gambl_metadata(seq_type_filter = "genome") %>% dplyr::filter(cohort == "DLBCL_cell_lines") %>% pull(sample_id)
  expect_false(all(get_coding_ssm(exclude_cohort = "DLBCL_cell_lines")[,"Tumor_Sample_Barcode"] %in% dlbcl_cell_lines_samples))
  
  #only return SSM for DLBCL_cell_line samples (cohort)
  expect_equal(get_coding_ssm(limit_cohort = "DLBCL_cell_lines"),
               get_coding_ssm(these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome") %>% 
                                dplyr::filter(cohort == "DLBCL_cell_lines")))
  
  #exclude DLBCL samples (pathology)
  dlbcl_samples = get_gambl_metadata(seq_type_filter = "genome") %>% dplyr::filter(pathology == "DLBCL") %>% pull(sample_id)
  dlbcl_ssm = get_coding_ssm(limit_pathology = "DLBCL") %>% pull(Tumor_Sample_Barcode)
  expect_true(all(dlbcl_ssm %in% dlbcl_samples))
  
  #use limit_samples to restrict the return to specific sample IDs
  expect_equal(get_coding_ssm(limit_samples = c("DOHH-2")), 
               get_coding_ssm(these_samples_metadata = get_gambl_metadata("genome") %>% 
                                dplyr::filter(sample_id %in% c("DOHH-2")))) 
})


test_that("Force unmatched samples", {
  expect_false(isTRUE(all.equal(get_coding_ssm(force_unmatched_samples = "DOHH-2", 
                                               these_samples_metadata = get_gambl_metadata("genome") %>% 
                                                 dplyr::filter(cohort == "DLBCL_cell_lines")),
                                
                                get_coding_ssm(these_samples_metadata = get_gambl_metadata("genome") %>% 
                                                 dplyr::filter(cohort == "DLBCL_cell_lines")))))
})


test_that("Check if `(this_seq_type = capture)` is working as inteded, i.e only capture samples are returned", {
  #get samples for each seq_type
  cap_samples = unique(get_gambl_metadata(seq_type_filter = "capture") %>% 
                         pull(sample_id))
  
  gen_samples = unique(get_gambl_metadata(seq_type_filter = "genome") %>% 
                         pull(sample_id))
  
  #request capture samples for tested get function
  cap_ssm = unique(get_coding_ssm(this_seq_type = "capture") %>% 
                     pull(Tumor_Sample_Barcode))
  
  #are all the requested samples in fact capture samples?
  expect_true(all(cap_ssm %in% cap_samples))
  
  #are the same sample set found in the genome sample pool?
  expect_false(all(cap_ssm %in% gen_samples))
})


test_that("Test if we can request specific unix groups", {
  gambl_samples = unique(get_gambl_metadata() %>% 
                           dplyr::filter(unix_group == "gambl") %>%
                           pull(sample_id))
  
  icgc_samples = unique(get_gambl_metadata() %>% 
                          dplyr::filter(unix_group == "icgc_dart") %>%
                          pull(sample_id))
  
  all_samples = unique(get_gambl_metadata() %>% 
                         dplyr::filter(unix_group %in% c("gambl", "icgc_dart")) %>%
                         pull(sample_id))

  gambl_ssm = unique(get_coding_ssm(groups = "gambl") %>% 
                     pull(Tumor_Sample_Barcode))
  
  icgc_ssm = unique(get_coding_ssm(groups = "icgc_dart") %>% 
                       pull(Tumor_Sample_Barcode))
  
  all_ssm = unique(get_coding_ssm(groups = c("gambl", "icgc_dart")) %>% 
                       pull(Tumor_Sample_Barcode))
  
  expect_true(all(gambl_ssm %in% gambl_samples))
  expect_true(all(icgc_ssm %in% icgc_samples))
  expect_true(all(all_ssm %in% all_samples))
})


test_that("See if other engines can be used (they can not, currently only fread_maf is the only supported engine...)", {
  expect_error(get_coding_ssm(engine = "V8"))
})


test_that("Test if the limit_samples is indeed deprecated", {
  expect_error(get_coding_ssm(limit_samples = "DOHH-2"))
})


test_that("Can we use these_sample_ids to return variants for a specifc set of samples", {
  expect_equal(get_coding_ssm(these_sample_ids = "DOHH-2"), 
               get_coding_ssm(these_samples_metadata = get_gambl_metadata() %>% 
                                dplyr::filter(sample_id == "DOHH-2")))
  
  expect_true(any(grepl("DOHH-2", get_coding_ssm(these_sample_ids = "DOHH-2") %>% 
                          pull(Tumor_Sample_Barcode) %>% 
                          unique())))
})
