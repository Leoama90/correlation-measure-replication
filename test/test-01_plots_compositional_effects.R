# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- test the that the preprocessed file is in the correct position --------

test_that("The file is not in the folder", {
  
  expected_data_frame <- here("script", "compositional_effects", "compositional_effects_02.rds")
  
  found_data_frame <- list.files(
    path       = here(),
    pattern    = "^compositional_effects_02\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  expect_setequal(found_data_frame, expected_data_frame)
})

# -------- test correct column name and size --------

test_that("there is a problem with the file columns: either one is empty or doesn't exist!", {
  df <- readRDS(here("script", "compositional_effects", "compositional_effects_02.rds"))
  
  expect_true(all(c("d", "pielou", "ERR_L1", "ERR_CLR") %in% names(df)))
  expect_gt(nrow(df), 0)
})


# -------- check the data frame after the pipelines --------

test_that("there is a problem with the pipelined dataframe!", {
  df <- readRDS(here("script", "compositional_effects", "compositional_effects_02.rds"))
  expect_true(nrow(df) > 0)
})


# -------- test pielou values are in [0, 1] --------

test_that("pielou values are outside [0, 1]", {
  df <- readRDS(here("script", "compositional_effects", "compositional_effects_02.rds"))
  
  expect_true(all(df$pielou >= 0 & df$pielou <= 1))
})


# -------- test d values are positive integers --------

test_that("d values are not positive integers", {
  df <- readRDS(here("script", "compositional_effects", "compositional_effects_02.rds"))
  
  expect_true(all(df$d > 0))
  expect_true(all(df$d == floor(df$d)))
})


# -------- test no NA in critical columns --------

test_that("there are NA values in critical columns", {
  df <- readRDS(here("script", "compositional_effects", "compositional_effects_02.rds"))
  
  expect_false(anyNA(df$d))
  expect_false(anyNA(df$pielou))
  expect_false(anyNA(df$ERR_L1))
  expect_false(anyNA(df$ERR_CLR))
})


# -------- test ERR_L1 and ERR_CLR are non-negative --------

test_that("ERR_L1 or ERR_CLR contain negative values", {
  df <- readRDS(here("script", "compositional_effects", "compositional_effects_02.rds"))
  
  expect_true(all(df$ERR_L1  >= 0))
  expect_true(all(df$ERR_CLR >= 0))
})


# -------- build df_sort (needed for downstream tests) --------

df      <- readRDS(here("script", "compositional_effects", "compositional_effects_02.rds"))
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


# -------- test df_sort dimensions --------

test_that("df_sort does not have the expected grid dimensions", {
  expected_rows <- length(seq(5, 200, by = 5)) * length(seq(0.025, 0.975, by = 0.025))
  
  expect_equal(nrow(df_sort), expected_rows)
})


# -------- test uniqueness of (d, pielou_round) combinations --------

test_that("df_sort contains duplicate (d, pielou_round) combinations", {
  n_distinct_pairs <- df_sort %>%
    dplyr::distinct(d, pielou_round) %>%
    nrow()
  
  expect_equal(n_distinct_pairs, nrow(df_sort))
})


# -------- test pielou_round contains exactly the expected grid values --------

test_that("pielou_round does not match the expected grid values", {
  expected_grid <- seq(0.025, 0.975, by = 0.025)
  
  expect_setequal(unique(df_sort$pielou_round), expected_grid)
})


# -------- test pielou_error is non-negative --------

test_that("pielou_error contains negative values", {
  expect_true(all(df_sort$pielou_error >= 0))
})


# -------- test pielou_error_logical is logical with no NA --------

test_that("pielou_error_logical is not logical or contains NA", {
  expect_type(df_sort$pielou_error_logical, "logical")
  expect_false(anyNA(df_sort$pielou_error_logical))
})


# -------- test df_percentiles has one row per unique d --------

test_that("df_percentiles does not have one row per unique value of d", {
  expect_equal(nrow(df_percentiles), length(unique(df_sort$d)))
})


# -------- test df_control proportion is acceptable --------

test_that("too many rows failed the pielou quality control threshold", {
  prop_failed <- nrow(df_control) / nrow(df_sort)
  
  expect_lt(prop_failed, 0.10)
})