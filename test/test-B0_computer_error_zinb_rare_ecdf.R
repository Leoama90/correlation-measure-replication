# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# load the script under test
source(here("script", "sparsity_effects", "B0_computer_error_zinb_rare_ecdf.R"))


# ---- OTU filtering & rarefaction -------------------------------------------

test_that("Column number after OTU filtering is not correct", {
  # prevalence and median filtering must reduce the number of OTUs
  expect_true(ncol(otu_filt) < ncol(otu_69001_H))
})

test_that("OTU prevalence filter was not applied correctly", {
  # every remaining OTU must appear in at least 33% of samples
  prevalences <- colSums(otu_filt > 0) / nrow(otu_filt)
  expect_true(all(prevalences >= 0.33))
})

test_that("OTU median filter was not applied correctly", {
  # every remaining OTU must have median non-zero abundance >= 5
  medians <- apply(otu_filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})

test_that("Rarefaction did not equalize library sizes", {
  # all samples must be rescaled to the same total read count
  row_sums <- rowSums(otu_filt)
  expect_true(all(row_sums == row_sums[1]))
})


# ---- ZINB parameter fitting -------------------------------------------------

test_that("HMP2_params does not have one row per filtered OTU", {
  # fitdistr must return one parameter set per OTU
  expect_equal(nrow(HMP2_params), ncol(otu_filt))
})

test_that("ZINB fitted parameters are outside their valid ranges", {
  # munb and size must be strictly positive; pstr0 must be a valid probability
  expect_true(all(HMP2_params$munb  >  0))
  expect_true(all(HMP2_params$size  >  0))
  expect_true(all(HMP2_params$pstr0 >= 0 & HMP2_params$pstr0 < 1))
})


# ---- Quantile filtering of ZINB parameters ----------------------------------

test_that("HMP2_params_filt retains OTUs outside the [10%, 90%] parameter ranges", {
  # all three parameters must lie within the computed decile bounds
  expect_true(all(HMP2_params_filt$munb  >= HMP2_quantile_params["10%", "munb" ] &
                    HMP2_params_filt$munb  <= HMP2_quantile_params["90%", "munb" ]))
  expect_true(all(HMP2_params_filt$size  >= HMP2_quantile_params["10%", "size" ] &
                    HMP2_params_filt$size  <= HMP2_quantile_params["90%", "size" ]))
  expect_true(all(HMP2_params_filt$pstr0 >= HMP2_quantile_params["10%", "pstr0"] &
                    HMP2_params_filt$pstr0 <= HMP2_quantile_params["90%", "pstr0"]))
})


# ---- Parameter grid (params_set) --------------------------------------------

test_that("params_set does not cover all expected cor x pstr0 combinations", {
  # grid must span 19 correlation values × 20 zero-inflation values = 380 rows
  expect_equal(nrow(params_set), length(seq(-.9, .9, by = .1)) *
                 length(seq(0, .95, by = .05)))
})


# ---- Final accumulated result -----------------------------------------------

test_that("result has wrong number of rows after all iterations", {
  # one row per params_set combination per outer iteration
  expect_equal(nrow(result), nIteration * nrow(params_set))
})

test_that("result cor_NorTA_PCLR contains invalid Pearson correlations", {
  # CLR-based Pearson correlations must be bounded in [-1, 1] with no NAs
  expect_false(anyNA(result$cor_NorTA_PCLR))
  expect_true(all(result$cor_NorTA_PCLR >= -1 & result$cor_NorTA_PCLR <= 1))
})