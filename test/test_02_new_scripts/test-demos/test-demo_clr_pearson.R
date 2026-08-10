# test-demo_clr_pearson.R
#
# Purpose:
#   This script tests the fundamental features of demo_clr_pearson.R,
#   checking that the pipeline runs end to end, that its result has the
#   expected structure, that the three output files are actually saved
#   to disk, and that the reloaded objects match the in-memory results.
#
# Inputs:
#   - demo_clr_pearson.R (sourced below)
#   - otu_HMP2.rds must be present in the project, since demo_clr_pearson.R
#     depends on it to run
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   Unlike test-clr_pearson.R, this file does not re-test the internal
#   correctness of clr_on_data() (already covered there). It only tests
#   the demo script's own responsibilities: running the pipeline on real
#   data and saving/reloading its results.

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# override readline() so the interactive prompts inside filt_data() and
# pseudocount() (called indirectly by the sourced demo script) are answered
# automatically with a valid threshold (0.1), instead of looping forever
# when tests run non-interactively. Assigned explicitly into .GlobalEnv,
# since filt_data()/pseudocount() are defined there (via source()'s default
# envir), and a plain top-level assignment in a testthat file may instead
# land in a local test environment, which those functions would never see.
assign("readline", function(prompt = "") "0.1", envir = .GlobalEnv)

# source demo_clr_pearson.R: this executes the whole pipeline (loading the
# raw data, running clr_on_data(), saving the three .rds files, and
# reloading them), leaving all its objects available in this environment.
# capture.output() suppresses its cat() messages so the test output stays clean.
capture.output(
  source(
    list.files(
      path = here(),
      pattern = "^demo_clr_pearson\\.R$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
)

# restore the original readline(), so this override doesn't leak into
# other test files run later in the same session
rm(readline, envir = .GlobalEnv)


# -------- test 1: the raw OTU data was loaded --------

test_that("raw OTU data was loaded into the environment", {
  expect_true(exists("otu"))
})


# -------- test 2: the pipeline result has the expected structure --------

test_that("pipeline result has the expected elements", {
  expect_named(result, c("samp_filt", "y_clr", "cor_matrix"))
})


# -------- test 3: all three output files were saved to disk --------

test_that("all three output files were saved to the outputs/ folder", {
  for (name in names(result)) {
    expect_true(file.exists(here("outputs", paste0(name, ".rds"))))
  }
})


# -------- test 4: the reloaded objects match the in-memory results --------

test_that("reloaded _view objects match the in-memory pipeline result", {
  expect_equal(samp_filt_view, result$samp_filt)
  expect_equal(y_clr_view, result$y_clr)
  expect_equal(cor_matrix_view, result$cor_matrix)
})


# -------- test 5: the filtered table cannot have more OTUs than the raw input --------

test_that("filtering does not increase the number of OTUs", {
  expect_true(ncol(result$samp_filt) <= ncol(otu))
})