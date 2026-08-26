# test-data_process.R
#
# Purpose:
#   This script tests the main features of dataprocess.R, checking:
#       - the three expected raw .csv files are found in the project
#       - none of the three files are empty (size > 0)
#       - the set of .csv files found in 00_data/raw/ matches exactly
#         the expected file inventory, with no extra or missing files
#
# Inputs:
#   - meta_HMP2.csv, otu_HMP2_16S.csv, taxonomy_HMP2_16S.csv
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- locate the .csv original data --------

file_01 <- list.files(
  path = here(),
  pattern = "^meta_HMP2\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
file_02 <- list.files(
  path = here(),
  pattern = "^otu_HMP2_16S\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
file_03 <- list.files(
  path = here(),
  pattern = "^taxonomy_HMP2_16S\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)


# -------- file existence --------

test_that("all three raw .csv files are found", {
  expect_true(file.exists(file_01))
  expect_true(file.exists(file_02))
  expect_true(file.exists(file_03))
})


# -------- file size --------

test_that("none of the raw .csv files are empty", {
  expect_true(file.size(file_01) > 0)
  expect_true(file.size(file_02) > 0)
  expect_true(file.size(file_03) > 0)
})


# -------- exact file inventory --------

test_that(".csv files found match the expected file inventory exactly", {
  expected_files <- c("meta_HMP2.csv", "otu_HMP2_16S.csv", "taxonomy_HMP2_16S.csv")
  # search only inside data/raw/ to avoid picking up unrelated CSVs
  # elsewhere in the project; recursive = FALSE keeps the search flat
  found_files <- basename(list.files(
    path      = here("00_data", "raw"),
    pattern   = "\\.csv$",
    recursive = FALSE
  ))
  # checks that both vectors contain exactly the same elements, regardless of order
  expect_setequal(found_files, expected_files)
})