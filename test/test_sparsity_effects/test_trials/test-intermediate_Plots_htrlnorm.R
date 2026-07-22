# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)


# -------- test: data loading and dataframe structure --------

test_that("df is loaded correctly and contains the expected columns", {
  # locate the dataset .rds file
  path <- list.files(
    path       = here(),
    pattern    = "Sparsity_Effects_htrlnorm.rds",
    full.names = TRUE,
    recursive  = TRUE
  )
  # ensure exactly one matching file was found
  expect_length(path, 1)

  # load the dataset and derive three new columns
  df <- readRDS(path) %>%
    # absolute error between observed and NorTA-PCLR-estimated correlation
    mutate("ERR_CLR" = abs(cor_normal - cor_NorTA_PCLR)) %>%
    # arithmetic mean of the two sparsity parameters
    mutate("phi_mean" = 0.5 * (phi_1 + phi_2)) %>%
    # element-wise maximum of phi_1 and phi_2
    mutate("phi_max" = pmax(phi_1, phi_2)) %>%
    # element-wise minimum of phi_1 and phi_2
    mutate("phi_min" = pmin(phi_1, phi_2))

  # verify all required columns are present
  expect_true(all(c(
    "ERR_CLR", "phi_mean", "phi_max", "phi_min",
    "cor_input", "cor_normal", "cor_NorTA_PCLR",
    "phi_1", "phi_2"
  ) %in% names(df)))
  # ERR_CLR is a correlation error, so it must be at most 1
  expect_false(any(df$ERR_CLR > 1))
  # ERR_CLR is an absolute value, so it must be non-negative
  expect_true(all(df$ERR_CLR >= 0))
})


# -------- test: consistency of derived sparsity columns --------

test_that("phi_min <= phi_mean <= phi_max for every row", {
  # locate the dataset .rds file
  path <- list.files(
    path       = here(),
    pattern    = "Sparsity_Effects_htrlnorm.rds",
    full.names = TRUE,
    recursive  = TRUE
  )
  # load data and compute the three sparsity summary columns
  df <- readRDS(path) %>%
    mutate("phi_mean" = 0.5 * (phi_1 + phi_2)) %>%
    mutate("phi_max" = pmax(phi_1, phi_2)) %>%
    mutate("phi_min" = pmin(phi_1, phi_2))

  # the minimum must never exceed the mean
  expect_true(all(df$phi_min <= df$phi_mean))
  # the mean must never exceed the maximum
  expect_true(all(df$phi_mean <= df$phi_max))
})


# -------- test: color palette --------

test_that("myPalette returns the correct number of colors and valid hex values", {
  # build a continuous palette by reversing the 11-color Spectral RColorBrewer palette
  myPalette <-
    RColorBrewer::brewer.pal(n = 11, "Spectral") %>% rev() %>%
    # wrap it in a color-ramp function that can interpolate any number of colors
    grDevices::colorRampPalette()

  # request 12 interpolated colors from the palette
  colors_12 <- myPalette(12)
  # confirm exactly 12 colors were returned
  expect_length(colors_12, 12)
  # validate each color string by trying to parse it with col2rgb();
  # col2rgb() throws an error for any invalid color string
  expect_true(all(sapply(colors_12, function(x) {
    tryCatch(
      {
        col2rgb(x)
        TRUE
      },
      error = function(e) FALSE
    )
  })))
})


# -------- test: aggregation via reframe --------

test_that("reframe computes ERR_CLR_MEAN correctly for each group", {
  # locate the dataset file
  path <- list.files(
    path       = here(),
    pattern    = "Sparsity_Effects_htrlnorm.rds",
    full.names = TRUE,
    recursive  = TRUE
  )
  # load data and derive error and minimum-sparsity columns
  df <- readRDS(path) %>%
    mutate("ERR_CLR" = abs(cor_normal - cor_NorTA_PCLR)) %>%
    mutate("phi_min" = pmin(phi_1, phi_2))

  # compute mean error for each (cor_input, phi_min) group
  agg <- df %>%
    reframe(
      ERR_CLR_MEAN = mean(ERR_CLR),
      .by = c(cor_input, phi_min)
    )

  # mean absolute error must be non-negative
  expect_true(all(agg$ERR_CLR_MEAN >= 0))
  # mean absolute error must be at most 1 (bounded correlation scale)
  expect_true(all(agg$ERR_CLR_MEAN <= 1))
  # each (cor_input, phi_min) combination must appear exactly once
  expect_equal(nrow(agg), nrow(distinct(agg, cor_input, phi_min)))
})


# -------- test: predict_max returns output of the correct size --------

test_that("predict_max returns a numeric value for each input", {
  # locate the OTU count matrix file
  f1 <- list.files(
    path = here(),
    pattern = "otu_HMP2.rds",
    full.names = TRUE,
    recursive = TRUE
  )

  # locate the sample metadata file
  f2 <- list.files(
    path = here(),
    pattern = "meta_HMP2.rds",
    full.names = TRUE,
    recursive = TRUE
  )

  # abort if either file is missing or duplicated
  stopifnot(length(f1) == 1, length(f2) == 1)

  # load OTU counts and sample metadata
  otu <- readRDS(f1)
  meta <- readRDS(f2)

  # keep only samples from subject 69-001 in the Healthy group
  otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]
  # retain taxa present in at least 33% of samples
  otu_filt <- otu_69001_H[, colSums(otu_69001_H > 0) / nrow(otu_69001_H) >= 0.33]
  # further retain taxa with a positive-read median >= 5
  otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]

  # build a tibble with log-scale max and mean abundance per taxon
  df_mean_max <- tibble(
    # log-transformed maximum across samples
    "y" = apply(log(otu_filt + 1), 2, max),
    # log-transformed mean across samples
    "x" = apply(log(otu_filt + 1), 2, mean)
  )
  # fit a simple linear model: max ~ mean (log scale)
  model <- lm(y ~ x, data = df_mean_max)

  # function that predicts the maximum for new mean-abundance values,
  # adding observation-level noise drawn from the fitted residual distribution
  predict_max <- function(new_xs) {
    new_data <- data.frame(x = new_xs)
    # get point predictions from the linear model
    predicted_values <- predict(model, new_data, interval = "none")
    # estimate residual variance from the model
    residuals_variance <- sum(residuals(model)^2) / model$df.residual
    # number of observations used to fit the model
    n <- length(model$model$x)
    # mean of the predictor in the training data
    x_bar <- mean(model$model$x)
    # leverage of each new point (how far it is from the training mean)
    leverages <- 1 / n + ((new_data$x - x_bar)^2 / sum((model$model$x - x_bar)^2))
    # prediction standard error, including both model and residual uncertainty
    std_error_prediction <- sqrt(residuals_variance * (1 + leverages))
    # draw one random prediction per new observation
    rnorm(nrow(new_data), mean = predicted_values, sd = std_error_prediction)
  }

  # call predict_max with three test values
  result <- predict_max(c(1, 2, 3))
  # must return exactly one value per input
  expect_length(result, 3)
  # output must be numeric
  expect_true(is.numeric(result))
  # no NA values are allowed
  expect_false(any(is.na(result)))
})
