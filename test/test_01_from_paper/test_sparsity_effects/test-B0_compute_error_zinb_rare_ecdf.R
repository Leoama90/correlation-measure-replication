# test-B0_compute_error_zinb_rare_ecdf.R
#
# Purpose:
#   This script tests the fundamental features of
#   B0_computer_error_zinb_rare_ecdf.R, checking:
#       - prevalence/median OTU filtering reduces the OTU count, and
#         every surviving OTU satisfies both thresholds
#       - rarefaction equalizes library sizes across all samples
#       - HMP2_params has one fitted ZINB parameter row per OTU, with
#         parameters within their valid ranges (munb, size > 0;
#         0 <= pstr0 < 1)
#       - HMP2_params_filt only retains OTUs whose parameters fall
#         within the computed [10%, 90%] quantile range
#       - params_set covers every combination of the correlation and
#         zero-inflation (pstr0) grid
#       - the final accumulated result has one row per params_set
#         combination per outer iteration, with valid (non-NA,
#         bounded in [-1, 1]) CLR Pearson correlations
#
# Inputs:
#   - B0_computer_error_zinb_rare_ecdf.R (sourced below)
#   - otu_HMP2.rds, meta_HMP2.rds, required by the sourced script
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   This script DOES source B0_computer_error_zinb_rare_ecdf.R directly,
#   which re-runs the full simulation on every test run: 100 outer
#   iterations, each with a parallel foreach over 380 (cor x pstr0)
#   combinations (38,000 total simulated OTU pairs, 10^4 samples each).
#   This is very slow and resource-heavy for a test suite; a faster
#   alternative (following the pattern used elsewhere in this project)
#   would load the already-saved Sparsity_Effects_zinbin_rare_ecdf.rds
#   output instead, or a small subset of it, rather than re-running the
#   whole simulation from scratch.

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- load the script under test --------

source(
  list.files(
    path = here(),
    pattern = "^B0_computer_error_zinb_rare_ecdf\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- OTU filtering & rarefaction --------

test_that("column number after OTU filtering is correct", {
  # prevalence and median filtering must reduce the number of OTUs
  expect_true(ncol(otu_filt) < ncol(otu_69001_H))
})

test_that("OTU prevalence filter was applied correctly", {
  # every remaining OTU must appear in at least 33% of samples
  prevalences <- colSums(otu_filt > 0) / nrow(otu_filt)
  expect_true(all(prevalences >= 0.33))
})

test_that("OTU median filter was applied correctly", {
  # every remaining OTU must have median non-zero abundance >= 5
  medians <- apply(otu_filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})

test_that("rarefaction equalized library sizes", {
  # all samples must be rescaled to the same total read count
  row_sums <- rowSums(otu_filt)
  expect_true(all(row_sums == row_sums[1]))
})


# -------- ZINB parameter fitting --------

test_that("HMP2_params has one row per filtered OTU", {
  # fitdistr must return one parameter set per OTU
  expect_equal(nrow(HMP2_params), ncol(otu_filt))
})

test_that("ZINB fitted parameters are within their valid ranges", {
  # munb and size must be strictly positive; pstr0 must be a valid probability
  expect_true(all(HMP2_params$munb > 0))
  expect_true(all(HMP2_params$size > 0))
  expect_true(all(HMP2_params$pstr0 >= 0 & HMP2_params$pstr0 < 1))
})


# -------- quantile filtering of ZINB parameters --------

test_that("HMP2_params_filt only retains OTUs within the [10%, 90%] parameter ranges", {
  # all three parameters must lie within the computed decile bounds
  expect_true(all(HMP2_params_filt$munb >= HMP2_quantile_params["10%", "munb"] &
                    HMP2_params_filt$munb <= HMP2_quantile_params["90%", "munb"]))
  expect_true(all(HMP2_params_filt$size >= HMP2_quantile_params["10%", "size"] &
                    HMP2_params_filt$size <= HMP2_quantile_params["90%", "size"]))
  expect_true(all(HMP2_params_filt$pstr0 >= HMP2_quantile_params["10%", "pstr0"] &
                    HMP2_params_filt$pstr0 <= HMP2_quantile_params["90%", "pstr0"]))
})


# -------- parameter grid (params_set) --------

test_that("params_set covers all expected cor x pstr0 combinations", {
  # grid must span 19 correlation values x 20 zero-inflation values = 380 rows
  expect_equal(
    nrow(params_set),
    length(seq(-0.9, 0.9, by = 0.1)) * length(seq(0, 0.95, by = 0.05))
  )
})


# -------- final accumulated result --------

test_that("result has the correct number of rows after all iterations", {
  # one row per params_set combination per outer iteration
  expect_equal(nrow(result), nIteration * nrow(params_set))
})

test_that("result's cor_NorTA_PCLR contains only valid Pearson correlations", {
  # CLR-based Pearson correlations must be bounded in [-1, 1] with no NAs
  expect_false(anyNA(result$cor_NorTA_PCLR))
  expect_true(all(result$cor_NorTA_PCLR >= -1 & result$cor_NorTA_PCLR <= 1))
})