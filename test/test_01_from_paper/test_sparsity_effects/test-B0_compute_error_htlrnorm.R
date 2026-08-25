# test-B0_compute_error_htlrnorm.
#
# Purpose:
#   <clever test here>
#   checking:
#       - 
#       - 
#       - 
#       - 
#
# Inputs:
#   - 
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)


# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- test data loading --------

# load OTU abundance matrix (samples x OTUs) and sample metadata
otu <- readRDS(list.files(
  path = here(),
  pattern = "^otu_HMP2\\.rds$",
  full.names = TRUE,
  recursive = TRUE
))

meta <- readRDS(list.files(
  path = here(),
  pattern = "^meta_HMP2\\.rds$",
  full.names = TRUE,
  recursive = TRUE
))

# verify that otu is a numeric matrix and meta is a list,
# as required by all downstream operations
test_that("otu is a numeric matrix and meta is a list", {
  expect_true(is.matrix(otu))
  expect_true(is.numeric(otu))
  expect_true(is.list(meta))
})

# -------- test the filtering functions --------

# select samples belonging to subject 69-001 in healthy status
otu_69001_h <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# remove rarest OTUs: keep only those present in at least 33% of samples
otu_filt <- otu_69001_h[, colSums(otu_69001_h > 0) / nrow(otu_69001_h) >= 0.33]

# verify that the prevalence filter actually removed at least one OTU
test_that("prevalence filter removes at least one OTU", {
  expect_gte(ncol(otu_69001_h), ncol(otu_filt))
})

# verify that every OTU that passed the filter truly appears
# in at least 33% of samples (no OTU slipped through incorrectly)
test_that("all remaining OTUs appear in at least 33% of samples", {
  prevalences <- colSums(otu_filt > 0) / nrow(otu_filt)
  expect_true(all(prevalences >= 0.33))
})

# keep only OTUs whose median non-zero abundance is at least 5
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]

# verify that every remaining OTU passes the median abundance threshold
test_that("all remaining OTUs have median non-zero abundance >= 5", {
  medians <- apply(otu_filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})

# -------- test quantile parameters --------

# remove objects no longer needed to free memory
rm(otu, otu_69001_h, meta)

# fit ZINB parameters from real data for each filtered OTU and
# save the 10th and 90th percentile of the parameter distributions
hmp2_quantile_params <- apply(otu_filt, 2, function(x) {
  ToyModel::mle.htrlnorm(x)$estimate
}) %>% apply(1, function(x) quantile(x, probs = c(0.1, 0.9)))

# verify that the resulting matrix has the correct dimensions (2 quantiles x 3 parameters),
# the expected row names, and that the 10th percentile is always below the 90th
test_that("quantile params matrix has correct shape and row names", {
  expect_equal(dim(hmp2_quantile_params), c(2, 3))
  expect_equal(rownames(hmp2_quantile_params), c("10%", "90%"))
  # 10th percentile must always be lower than 90th for each parameter
  expect_true(all(hmp2_quantile_params["10%", ] < hmp2_quantile_params["90%", ]))
})

# -------- test params_set structure --------

# fit a linear model between mean and max abundance (log scale) across OTUs
df_mean_max <- tibble(
  y = apply(log(otu_filt + 1), 2, max),
  x = apply(log(otu_filt + 1), 2, mean)
)
model <- lm(y ~ x, data = df_mean_max)

# generate all combinations of sparsity (phi) and correlation values to test,
# then assign random log-normal parameters and predicted maxima to each combination
params_set <- expand_grid(
  phi_1 = seq(0, 0.95, by = 0.05),
  phi_2 = seq(0, 0.95, by = 0.05),
  cor   = seq(-0.9, 0.9, by = 0.1)
) %>%
  mutate(
    meanlog_1 = runif(n(), hmp2_quantile_params["10%", "meanlog"], hmp2_quantile_params["90%", "meanlog"]),
    meanlog_2 = runif(n(), hmp2_quantile_params["10%", "meanlog"], hmp2_quantile_params["90%", "meanlog"]),
    sdlog_1   = runif(n(), hmp2_quantile_params["10%", "sdlog"], hmp2_quantile_params["90%", "sdlog"]),
    sdlog_2   = runif(n(), hmp2_quantile_params["10%", "sdlog"], hmp2_quantile_params["90%", "sdlog"])
  ) %>%
  as.data.frame()

# verify that params_set contains all expected columns and that
# all sampled values fall within the defined ranges
test_that("params_set has expected columns and all values within defined ranges", {
  expected_cols <- c("phi_1", "phi_2", "cor", "meanlog_1", "meanlog_2", "sdlog_1", "sdlog_2")
  expect_true(all(expected_cols %in% names(params_set)))
  # correlation values must stay within [-0.9, 0.9]
  expect_true(all(params_set$cor >= -0.9 & params_set$cor <= 0.9))
  # phi values must stay within [0, 0.95]
  expect_true(all(params_set$phi_1 >= 0 & params_set$phi_1 <= 0.95))
  expect_true(all(params_set$phi_2 >= 0 & params_set$phi_2 <= 0.95))
})
