# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- Setup: load data once, shared across all tests --------

rds_path <- list.files(
  path       = here(),
  pattern    = "^compositional_effects_02\\.rds$",
  recursive  = TRUE,
  full.names = TRUE
)

stopifnot("compositional_effects_02.rds not found or duplicated" = length(rds_path) == 1)

df       <- readRDS(rds_path)

df_sort <- tibble::tibble()

for (di in seq(5, 200, by = 5)) {
  for (ei in seq(0.025, 0.975, by = 0.025)) {
    sub_df  <- df %>% dplyr::filter(d == di)
    row_i   <- sub_df[which.min(abs(sub_df$pielou - ei)), ]
    row_i   <- c(row_i, "pielou_round" = ei)
    df_sort <- dplyr::bind_rows(df_sort, row_i)
  }
}

df_sort <- df_sort %>%
  dplyr::mutate(pielou_error         = abs(pielou_round - pielou)) %>%
  dplyr::mutate(pielou_error_logical = pielou_error < 0.005, .after = pielou_error) %>%
  dplyr::mutate(log_err_l1           = log10(ERR_L1)) %>%
  dplyr::mutate(log_err_l1           = dplyr::if_else(log_err_l1  < -2, -2, log_err_l1)) %>%
  dplyr::mutate(log_err_clr          = log10(ERR_CLR)) %>%
  dplyr::mutate(log_err_clr          = dplyr::if_else(log_err_clr < -2, -2, log_err_clr))

df_percentiles <- df_sort %>%
  dplyr::group_by(d) %>%
  dplyr::summarise(
    err_clr_mean = mean(ERR_CLR),
    err_clr_p10  = quantile(ERR_CLR, 0.10),
    err_clr_p90  = quantile(ERR_CLR, 0.90)
  )

df_control <- df_sort %>% dplyr::filter(pielou_error_logical == FALSE)


# -------- test that the preprocessed file is in the expected location --------

test_that("compositional_effects_02.rds must be at the expected project path", {

  found_data_frame <- list.files(
    path       = here(),
    pattern    = "^compositional_effects_02\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  )

  expect_equal(found_data_frame, rds_path)
})


# -------- test correct column names and non-empty rows --------

test_that("df must contain columns d, pielou, ERR_L1, ERR_CLR and at least one row", {

  expect_true(all(c("d", "pielou", "ERR_L1", "ERR_CLR") %in% names(df)))
  expect_gt(nrow(df), 0)
})


# -------- test pielou values are in [0, 1] --------

test_that("pielou values must lie within [0, 1]", {

  expect_true(all(df$pielou >= 0 & df$pielou <= 1))
})


# -------- test d values are positive integers --------

test_that("d must contain only positive integers", {

  expect_true(all(df$d > 0))
  expect_true(all(df$d == floor(df$d)))
})


# -------- test no NA in critical columns --------

test_that("critical columns (d, pielou, ERR_L1, ERR_CLR) must contain no NA values", {

  expect_false(anyNA(df$d))
  expect_false(anyNA(df$pielou))
  expect_false(anyNA(df$ERR_L1))
  expect_false(anyNA(df$ERR_CLR))
})


# -------- test ERR_L1 and ERR_CLR are non-negative --------

test_that("ERR_L1 and ERR_CLR must be non-negative", {

  expect_true(all(df$ERR_L1  >= 0))
  expect_true(all(df$ERR_CLR >= 0))
})


# -------- test df_sort dimensions --------

test_that("df_sort must have one row per (d, pielou_round) grid combination", {

  expected_rows <- length(seq(5, 200, by = 5)) * length(seq(0.025, 0.975, by = 0.025))

  expect_equal(nrow(df_sort), expected_rows)
})


# -------- test uniqueness of (d, pielou_round) combinations --------

test_that("each (d, pielou_round) pair in df_sort must be unique", {

  n_distinct_pairs <- df_sort %>%
    dplyr::distinct(d, pielou_round) %>%
    nrow()

  expect_equal(n_distinct_pairs, nrow(df_sort))
})


# -------- test pielou_round contains exactly the expected grid values --------

test_that("pielou_round must contain exactly the values of the expected grid seq(0.025, 0.975, by=0.025)", {

  expected_grid <- seq(0.025, 0.975, by = 0.025)

  expect_setequal(unique(df_sort$pielou_round), expected_grid)
})


# -------- test pielou_error is non-negative --------

test_that("pielou_error must be non-negative (it is an absolute difference)", {

  expect_true(all(df_sort$pielou_error >= 0))
})


# -------- test pielou_error_logical is logical with no NA --------

test_that("pielou_error_logical must be of type logical and contain no NA values", {

  expect_type(df_sort$pielou_error_logical, "logical")
  expect_false(anyNA(df_sort$pielou_error_logical))
})


# -------- test df_percentiles has one row per unique d --------

test_that("df_percentiles must have exactly one summary row per unique value of d", {

  expect_equal(nrow(df_percentiles), length(unique(df_sort$d)))
})


# -------- test df_control proportion is acceptable --------

test_that("rows failing the pielou quality-control threshold must be fewer than 10%", {

  prop_failed <- nrow(df_control) / nrow(df_sort)

  expect_lt(prop_failed, 0.10)
})