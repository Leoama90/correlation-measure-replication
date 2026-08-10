# test-filt_data.R
#
# Purpose:
#   This script tests the fundamental features of filt_data(), checking
#   that OTUs are correctly kept or dropped based on the prevalence and
#   median-abundance thresholds, that the number of samples (rows) is
#   left unchanged, and that the optional taxa_filt output is included
#   only when a "taxa" table exists in the calling environment.
#
# Inputs:
#   - filt_data.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   filt_data() now accepts prevalence_threshold and abundance_threshold
#   as explicit arguments; passing both avoids the interactive readline()
#   prompt entirely, so the actual function can be tested directly,
#   without needing a non-interactive duplicate of its body.

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# bring filt_data() (and datasum(), sourced inside it) into scope
source(here("02_new_scripts", "filt_data.R"))

# suppress cat() output (from datasum(), called inside filt_data())
# so tests stay clean, keep the return value
quiet_filt_data <- function(x, prevalence_threshold, abundance_threshold = 5) {
  capture.output(
    result <- filt_data(
      x,
      prevalence_threshold = prevalence_threshold,
      abundance_threshold = abundance_threshold
    )
  )
  result
}


# -------- shared test data --------

# a small OTU table (5 samples x 3 OTUs), built so each OTU is dropped
# or kept for a distinct, known reason at the 0.1 threshold used below:
#   - OTU1: present in every sample, decent counts -> should be kept
#   - OTU2: present in no sample (all zero) -> dropped by prevalence
#   - OTU3: present in every sample, but very low counts -> dropped by
#     the median (>= 5 reads) filter
otu_test <- matrix(
  c(
    10, 0, 1,
    12, 0, 2,
    11, 0, 1,
     9, 0, 1,
    13, 0, 2
  ),
  nrow = 5, byrow = TRUE,
  dimnames = list(NULL, c("OTU1", "OTU2", "OTU3"))
)


# -------- test 1: OTUs are correctly kept/dropped by the two filters --------

test_that("filtering keeps only the OTU passing both prevalence and median filters", {
  # make sure no leftover "taxa" object affects this test
  if (exists("taxa")) rm(taxa, envir = .GlobalEnv)
  result <- quiet_filt_data(otu_test, prevalence_threshold = 0.1)
  # only OTU1 satisfies both the prevalence and the median-abundance filter
  expect_equal(colnames(result$samp_filt), "OTU1")
})


# -------- test 2: filtering does not change the number of samples --------

test_that("filtering leaves the number of samples (rows) unchanged", {
  if (exists("taxa")) rm(taxa, envir = .GlobalEnv)
  result <- quiet_filt_data(otu_test, prevalence_threshold = 0.1)
  # filtering only removes columns (OTUs), never rows (samples)
  expect_equal(nrow(result$samp_filt), nrow(otu_test))
})


# -------- test 3: taxa_filt is absent when no taxa table exists --------

test_that("result has no taxa_filt element when 'taxa' does not exist", {
  if (exists("taxa")) rm(taxa, envir = .GlobalEnv)
  result <- quiet_filt_data(otu_test, prevalence_threshold = 0.1)
  # the taxa branch is only triggered if a "taxa" object is found
  expect_null(result$taxa_filt)
})


# -------- test 4: taxa_filt is included and aligned when taxa exists --------

test_that("result includes taxa_filt, aligned to the surviving OTUs, when 'taxa' exists", {
  # define a taxa table in the global environment, indexed by OTU name
  taxa <<- data.frame(
    genus = c("GenusA", "GenusB", "GenusC"),
    row.names = c("OTU1", "OTU2", "OTU3")
  )
  result <- quiet_filt_data(otu_test, prevalence_threshold = 0.1)
  # taxa_filt should contain only the row matching the surviving OTU (OTU1)
  expect_equal(rownames(result$taxa_filt), "OTU1")
  # clean up, so later tests are not affected by this global object
  rm(taxa, envir = .GlobalEnv)
})


# -------- test 5: abundance_threshold parameter is respected --------

test_that("a lower abundance_threshold rescues an OTU otherwise dropped by the median filter", {
  if (exists("taxa")) rm(taxa, envir = .GlobalEnv)
  # OTU3 has non-zero values of 1 and 2 (median well below the default
  # threshold of 5); lowering abundance_threshold to 1 should let it pass
  result <- quiet_filt_data(otu_test, prevalence_threshold = 0.1, abundance_threshold = 1)
  expect_true("OTU3" %in% colnames(result$samp_filt))
})
