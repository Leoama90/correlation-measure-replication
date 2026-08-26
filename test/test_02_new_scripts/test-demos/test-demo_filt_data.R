# test-demo_filt_data.R
#
# Purpose:
#   Tests that running demo_filt_data.R produces the expected objects
#   and results: a correctly-shaped example OTU table, a filtered
#   result consistent with the known lambda values used to generate
#   it, and a correctly aligned taxa_filt when a taxa table is present.
#   This does not re-test filt_data() itself (see test-filt_data.R),
#   only that the demo script uses it correctly.
#
# Inputs:
#   - demo_filt_data.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# run the demo script with its cat()/print() output suppressed, keeping
# the objects it creates (otu_demo, filt_result, filt_result_with_taxa)
# available in this environment for the tests below
capture.output(
  source(
    list.files(
      path = here(),
      pattern = "^demo_filt_data\\.R$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
)

# -------- test 1: the example OTU table has the expected shape --------

test_that("otu_demo has the expected dimensions and OTU names", {
  expect_equal(dim(otu_demo), c(10, 5))
  expect_equal(colnames(otu_demo), paste0("OTU", 1:5))
})


# -------- test 2: filtering drops the low-prevalence/low-abundance OTUs --------

test_that("filt_result keeps only the OTUs expected from the known lambda values", {
  # OTU1, OTU3, OTU4 have high enough lambda to pass both filters;
  # OTU2 (lambda 0.5) is dropped by prevalence, OTU5 (lambda 1) by
  # the abundance/median filter
  expect_true(all(c("OTU1", "OTU3", "OTU4") %in% colnames(filt_result$samp_filt)))
  expect_false("OTU2" %in% colnames(filt_result$samp_filt))
  expect_false("OTU5" %in% colnames(filt_result$samp_filt))
})


# -------- test 3: taxa_filt is aligned to the surviving OTUs --------

test_that("filt_result_with_taxa$taxa_filt matches the OTUs kept after filtering", {
  expect_equal(
    rownames(filt_result_with_taxa$taxa_filt),
    colnames(filt_result_with_taxa$samp_filt)
  )
})
