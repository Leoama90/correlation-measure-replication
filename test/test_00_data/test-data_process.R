# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- file existence --------

test_that("file not found, where did I put it...?", {
  # give the file position, here() searches it, then store in vector file_01
  file_01 <- here("data", "raw", "meta_HMP2.csv")
  # logic, if file.exists = TRUE, test pass
  expect_true(file.exists(file_01))

  # give the file position, here() searches it, then store in vector file_02
  file_02 <- here("data", "raw", "otu_HMP2_16S.csv")
  # logic, if file.exists = TRUE, test pass
  expect_true(file.exists(file_02))

  # give the file position, here() searches it, then store in vector file_03
  file_03 <- here("data", "raw", "taxonomy_HMP2_16S.csv")
  # logic, if file.exists = TRUE, test pass
  expect_true(file.exists(file_03))
})


# -------- file size --------

test_that("file size = 0, maybe it's empty and that's a problem!", {
  # test that each file is not empty (size > 0)
  file_01 <- here("data", "raw", "meta_HMP2.csv")
  expect_true(file.size(file_01) > 0)

  file_02 <- here("data", "raw", "otu_HMP2_16S.csv")
  expect_true(file.size(file_02) > 0)

  file_03 <- here("data", "raw", "taxonomy_HMP2_16S.csv")
  expect_true(file.size(file_03) > 0)
})


# -------- exact file inventory --------

test_that(".csv files found are NOT the expected ones", {
  expected_files <- c("meta_HMP2.csv", "otu_HMP2_16S.csv", "taxonomy_HMP2_16S.csv")

  # search only inside data/raw/ to avoid picking up unrelated CSVs
  # elsewhere in the project; recursive = FALSE keeps the search flat
  found_files <- basename(list.files(
    path      = here("data", "raw"),
    pattern   = "\\.csv$",
    recursive = FALSE
  ))

  # checks that both vectors contain exactly the same elements, regardless of order
  expect_setequal(found_files, expected_files)
})
