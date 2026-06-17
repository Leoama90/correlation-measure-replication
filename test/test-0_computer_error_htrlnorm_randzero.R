# doSNOW: parallel backend for foreach, supports progress bars via snow clusters
# https://cran.r-project.org/web/packages/doSNOW/index.html
library(doSNOW)

# foreach: parallel foreach loops
# https://cran.r-project.org/package=foreach
library(foreach)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# progress: displays text progress bars for long-running loops
# https://cran.r-project.org/web/packages/progress/index.html
library(progress)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# ToyModel: lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)


# -------- Important Consideration! --------

# for this test, I will load some chunk of code taken from the original one,
# because loading the whole script would a too long execution time


# -------- 1. data filtering code --------

# load the first chunk of code from the script 0_computer_error_htrlnorm_randzero.R
# that will allow the creation of the test

# Load OTU abundance matrix (samples x OTUs)
otu  <- readRDS(here("data", "otu_HMP2.rds" ))
# Load sample metadata
meta <- readRDS(here("data", "meta_HMP2.rds"))

# Subset samples belonging to subject 69-001 in the healthy state
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# Remove OTUs present in fewer than 33% of samples (low prevalence)
otu_filt <- otu_69001_H[, colSums(otu_69001_H > 0)/nrow(otu_69001_H) >= 0.33]
# Further remove OTUs whose median non-zero count is below 5 (too rare)
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]


# -------- test the 1. chunk --------

test_that("The row number is not correct", {
  # get otu row number
  otu_rown    <- nrow(otu)
  # get otu_69001_H row number
  otu_69_rown <- nrow(otu_69001_H)
  # expect less rows after filtering for subject and health status
  expect_true(otu_rown > otu_69_rown)
})

test_that("The column number is not correct", {
  # expect less OTUs after prevalence and median filtering
  expect_true(ncol(otu_filt) < ncol(otu_69001_H))
})


# -------- 2. linear model --------

# Build a data frame pairing each OTU's log-mean with its log-maximum,
# used to learn a realistic mean–max relationship from the real data
df_mean_max <- tibble(
  "y" = apply(log(otu_filt + 1), 2, max ),
  "x" = apply(log(otu_filt + 1), 2, mean)
)

# Fit a simple linear model: log(max) ~ log(mean)
model <- lm(y~x, data = df_mean_max)


# -------- test the 2. chunk --------

# create a dummy matrix to use it to test the linear model

make_otu <- function(nrow = 10, ncol = 5, seed = 42) {
  set.seed(seed)
  m <- matrix(rpois(nrow * ncol, lambda = 10), 
              nrow = nrow, 
              ncol = ncol)
  colnames(m) <- paste0("OTU", seq_len(ncol))
  m
}

# Each OTU becomes one row, so row count must equal the number of OTU columns
test_that("df_mean_max has as many rows as columns in otu_filt", {
  otu_filt <- make_otu(ncol = 7)
  df_mean_max <- tibble(
    y = apply(log(otu_filt + 1), 2, max),
    x = apply(log(otu_filt + 1), 2, mean)
  )
  expect_equal(nrow(df_mean_max), ncol(otu_filt))
})

# The maximum of a set of values is always >= its mean, so y >= x must hold
test_that("y (log-max) is always >= x (log-mean) for every OTU", {
  otu_filt <- make_otu()
  df_mean_max <- tibble(
    y = apply(log(otu_filt + 1), 2, max),
    x = apply(log(otu_filt + 1), 2, mean)
  )
  expect_true(all(df_mean_max$y >= df_mean_max$x))
})

# log(count + 1) with count >= 0 guarantees non-negative transformed values
test_that("x and y are all >= 0 due to the log(count + 1) shift", {
  otu_filt <- make_otu()
  df_mean_max <- tibble(
    y = apply(log(otu_filt + 1), 2, max),
    x = apply(log(otu_filt + 1), 2, mean)
  )
  expect_true(all(df_mean_max$x >= 0))
  expect_true(all(df_mean_max$y >= 0))
})

# The stored formula must be exactly y ~ x, not a rearranged or expanded form
test_that("the model formula is y ~ x", {
  otu_filt <- make_otu()
  df_mean_max <- tibble(
    y = apply(log(otu_filt + 1), 2, max),
    x = apply(log(otu_filt + 1), 2, mean)
  )
  model <- lm(y ~ x, data = df_mean_max)
  expect_equal(deparse(formula(model)), "y ~ x")
})

# predict() must return one fitted value per row in the input data frame
test_that("predict() returns as many values as rows in df_mean_max", {
  otu_filt <- make_otu()
  df_mean_max <- tibble(
    y = apply(log(otu_filt + 1), 2, max),
    x = apply(log(otu_filt + 1), 2, mean)
  )
  model <- lm(y ~ x, data = df_mean_max)
  preds <- predict(model, newdata = df_mean_max)
  expect_length(preds, nrow(df_mean_max))
})


# -------- test for both inner and outer loop --------

# Fixed HMP2 quantile parameter boundaries used throughout the tests
HMP2_quantile_params <- matrix(
  c(1.0, 0.3, 0.1,
    3.0, 1.2, 0.8),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("10%", "90%"), c("meanlog", "sdlog", "phi"))
)

n_background_otus <- 20
n_samples         <- 100


# -------- test for the outer loop --------

# params_random_HMP2 must have one row per OTU and values within HMP2 bounds
test_that("params_random_HMP2 has correct shape and values within HMP2 bounds", {
  set.seed(1)
  params_random_HMP2 <- data.frame(
    meanlog = runif(n_background_otus,
                    HMP2_quantile_params["10%", "meanlog"],
                    HMP2_quantile_params["90%", "meanlog"]),
    sdlog   = runif(n_background_otus,
                    HMP2_quantile_params["10%", "sdlog"],
                    HMP2_quantile_params["90%", "sdlog"]),
    phi     = runif(n_background_otus,
                    HMP2_quantile_params["10%", "phi"],
                    HMP2_quantile_params["90%", "phi"])
  )
  expect_equal(nrow(params_random_HMP2), n_background_otus)
  expect_named(params_random_HMP2, c("meanlog", "sdlog", "phi"))
  expect_true(all(params_random_HMP2$meanlog >= HMP2_quantile_params["10%", "meanlog"]))
  expect_true(all(params_random_HMP2$meanlog <= HMP2_quantile_params["90%", "meanlog"]))
  expect_true(all(params_random_HMP2$phi     >= HMP2_quantile_params["10%", "phi"]))
  expect_true(all(params_random_HMP2$phi     <= HMP2_quantile_params["90%", "phi"]))
})

# params_set must be a full factorial grid with the correct row count and columns
test_that("params_set is a full factorial grid with all required columns", {
  phi_seq <- seq(0, 0.95, by = 0.05)  # 20 levels
  cor_seq <- seq(-0.9, 0.9, by = 0.1) # 19 levels
  set.seed(1)
  params_set <- expand_grid(
    phi_1 = phi_seq,
    phi_2 = phi_seq,
    cor   = cor_seq
  ) %>%
    mutate(meanlog_1 = runif(n(), 1  ,   3),
           meanlog_2 = runif(n(), 1  ,   3),
           sdlog_1   = runif(n(), 0.3, 1.2),
           sdlog_2   = runif(n(), 0.3, 1.2))
  # Row count must equal 20 × 20 × 19 = 7,600
  expect_equal(nrow(params_set), length(phi_seq)^2 * length(cor_seq))
  expect_true(all(c("phi_1", "phi_2", "cor",
                    "meanlog_1", "meanlog_2",
                    "sdlog_1",   "sdlog_2") %in% names(params_set)))
})


# -------- test for the inner loop --------

# After grafting two OTUs into the background matrix, only the target columns
# must change; all other columns must remain identical to the original
test_that("replacing target columns leaves all other columns unchanged", {
  set.seed(42)
  background <- matrix(rpois(n_samples * n_background_otus, 5),
                       nrow = n_samples, ncol = n_background_otus)
  new_col_a  <- rpois(n_samples, lambda = 8)
  new_col_b  <- rpois(n_samples, lambda = 3)
  modified        <- background
  modified[, 5]  <- new_col_a
  modified[, 15] <- new_col_b
  expect_equal(modified[, 5],  new_col_a)
  expect_equal(modified[, 15], new_col_b)
  # Every column other than 5 and 15 must be untouched
  expect_equal(modified[, -c(5, 15)], background[, -c(5, 15)])
})

# After pseudo-value imputation no exact zeros must remain, and non-zero
# entries must be unchanged
test_that("pseudo-value imputation removes all zeros without altering non-zero entries", {
  set.seed(42)
  mat <- matrix(c(0, 1, 2, 0, 5, 0), nrow = 2)
  rand_pseudo <- matrix(runif(length(mat), min = 0.065, max = 0.65),
                        nrow = nrow(mat), ncol = ncol(mat))
  mat_imputed <- mat + (mat == 0) * rand_pseudo
  expect_true(all(mat_imputed > 0))
  expect_equal(mat_imputed[mat != 0], mat[mat != 0])
})

# The CLR-based correlation matrix must be symmetric, have ones on the diagonal,
# and all off-diagonal values within [-1, 1]
test_that("CLR correlation matrix is valid (symmetric, diagonal = 1, values in [-1,1])", {
  set.seed(42)
  # Build a positive matrix (no zeros) so log is always defined
  mat     <- matrix(rpois(n_samples * n_background_otus, 10) + 1,
                    nrow = n_samples, ncol = n_background_otus)
  log_mat <- log(mat)
  clr_mat <- log_mat - rowMeans(log_mat)
  cor_mat <- cor(clr_mat)
  expect_equal(dim(cor_mat), c(n_background_otus, n_background_otus))
  expect_equal(diag(cor_mat), rep(1, n_background_otus))
  expect_true(isSymmetric(cor_mat))
  expect_true(all(cor_mat >= -1 & cor_mat <= 1))
})

# The one-row result data frame must contain all 12 documented columns and
# the correlation estimate must be a finite scalar within [-1, 1]
test_that("inner loop result row has all required columns and a valid correlation", {
  set.seed(42)
  n_col   <- 20
  mat     <- matrix(rpois(n_samples * n_col, 10) + 1,
                    nrow = n_samples, ncol = n_col)
  log_mat  <- log(mat)
  clr_mat  <- log_mat - rowMeans(log_mat)
  cor_pclr <- cor(clr_mat)
  row_result <- data.frame(
    iteration      = 1L,
    cor_input      = 0.5,
    cor_normal     = 0.48,
    meanlog_1      = 2.1,
    meanlog_2      = 1.9,
    sdlog_1        = 0.6,
    sdlog_2        = 0.7,
    b_1            = 4.0,
    b_2            = 3.8,
    phi_1          = 0.2,
    phi_2          = 0.3,
    cor_NorTA_PCLR = cor_pclr[5, 15]
  )
  expected_cols <- c("iteration", "cor_input", "cor_normal",
                     "meanlog_1", "meanlog_2", "sdlog_1", "sdlog_2",
                     "b_1", "b_2", "phi_1", "phi_2", "cor_NorTA_PCLR")
  expect_named(row_result, expected_cols)
  expect_length(row_result$cor_NorTA_PCLR, 1)
  expect_true(is.finite(row_result$cor_NorTA_PCLR))
  expect_true(abs(row_result$cor_NorTA_PCLR) <= 1)
})

# After N iterations the accumulated data frame must have N × K rows and no NAs
test_that("accumulated result has correct row count and no NAs after all iterations", {
  set.seed(42)
  n_iter   <- 3
  n_params <- 5
  result   <- data.frame()
  for (iter in seq_len(n_iter)) {
    dummy_df <- data.frame(
      iteration      = iter,
      cor_input      = runif(n_params, -0.9,  0.9 ),
      cor_normal     = runif(n_params, -0.9,  0.9 ),
      meanlog_1      = runif(n_params,  1  ,  3   ),
      meanlog_2      = runif(n_params,  1  ,  3   ),
      sdlog_1        = runif(n_params,  0.3,  1.2 ),
      sdlog_2        = runif(n_params,  0.3,  1.2 ),
      b_1            = runif(n_params,  2  ,  5   ),
      b_2            = runif(n_params,  2  ,  5   ),
      phi_1          = runif(n_params,  0  ,  0.95),
      phi_2          = runif(n_params,  0  ,  0.95),
      cor_NorTA_PCLR = runif(n_params, -1  ,  1   )
    )
    result <- bind_rows(result, dummy_df)
  }
  expect_equal(nrow(result), n_iter * n_params)
  expect_false(anyNA(result))
})


