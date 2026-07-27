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
#   filt_data() uses a readline() command which goes into an infinite loop
#   once the script is sourced in the test. To avoid this, I duplicated the
#   body of the function in the test without the problematic command.
#   the trade-off is that this test must be changed accordingly to the original
#   script, if the script changes in the future.

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load filt_data.R only for reference/documentation purposes; the
# actual function used in the tests below is the non-interactive
# duplicate defined right after
source(here("new_scripts", "filt_data.R"))

# -------- non-interactive duplicate of filt_data(), for testing --------
# same logic as filt_data.R, but prevalence_threshold is a direct
# argument instead of being requested via readline()
filt_data_for_test <- function(x, prevalence_threshold) {
  # show the "before" summary; datasum() prints its own stats and
  # invisibly returns them, so no extra cat() is needed around it
  datasum(x)
  # filter the OTU table by the given prevalence threshold
  # drop = FALSE keeps the result as a data.frame/matrix even if only
  # one OTU survives the filter
  samp_filt <- x[, colSums(x > 0) / nrow(x) >= prevalence_threshold, drop = FALSE]
  # -------- median filter --------
  # number of OTUs still in the table after the prevalence filter
  n_otu <- ncol(samp_filt)
  # logical vector marking which OTUs pass the median (>= 5 reads) filter
  keep <- logical(n_otu)
  for (i in seq_len(n_otu)) {
    col <- samp_filt[, i]
    keep[i] <- median(col[col > 0]) >= 5
  }
  # keep only the OTUs whose median non-zero abundance is >= 5 reads
  samp_filt <- samp_filt[, keep, drop = FALSE]
  # show the "after" summary
  datasum(samp_filt)
  # if a taxa table exists in the calling environment, align it to
  # the OTUs that survived filtering; otherwise skip this step
  result <- list(samp_filt = samp_filt)
  if (exists("taxa")) {
    taxa_filt <- taxa[colnames(samp_filt), , drop = FALSE]
    result$taxa_filt <- taxa_filt
  }
  # return result as an invisible value
  invisible(result)
}

# suppress cat() output (from datasum(), called inside the duplicate)
# so tests stay clean, keep the return value
quiet_filt_data <- function(x, prevalence_threshold) {
  capture.output(result <- filt_data_for_test(x, prevalence_threshold))
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
  c(10, 0, 1,
    12, 0, 2,
    11, 0, 1,
    9, 0, 1,
    13, 0, 2),
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