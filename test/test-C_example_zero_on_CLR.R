# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# Load CLR function
source(here("script", "method_comparison", "CLR.R"))

# Set seed for reproducibility
set.seed(42)

# Define number of samples and OTUs to create a dummy matrix
n_samples <- 20
n_otus    <- 10

# new function that builds a random sparse OTU count matrix
make_otu <- function() {
  # Fill an n_samples x n_otus matrix with random integers in [0, 100]
  m <- matrix(
    sample(0:100, n_samples * n_otus, replace = TRUE),
    nrow = n_samples,
    # Assign row and column names for samples and OTUs
    dimnames = list(
      paste0("Sample_", seq_len(n_samples)),
      paste0("OTU_",    seq_len(n_otus))
    )
  )
  # Set 40% of entries to zero to simulate sparsity
  m[sample(length(m), size = floor(length(m) * 0.4))] <- 0
  m
}

# Build the dummy OTU matrix to use in all other tests
otu_raw <- make_otu()


# Check CLR output matrix has same number of rows and columns as the input
test_that("CLR output has same dimensions as input", {
  # Add a pseudocount to avoid log(0) in the CLR transformation
  otu_pos <- otu_raw + 0.5
  # Apply the CLR transformation
  otu_clr <- CLR(otu_pos)
  
  # Verify that input and output dimensions match exactly
  expect_equal(dim(otu_clr), dim(otu_pos))
})

# Check that each sample's CLR values sum to zero, as required by the CLR definition
test_that("CLR rows sum to (approximately) zero", {
  otu_pos <- otu_raw + 0.5
  otu_clr <- CLR(otu_pos)
  
  # Compute the row sum of CLR values
  row_sums <- rowSums(otu_clr)
  # Accept numerical near-zero (floating point tolerance 1e-10)
  expect_true(all(abs(row_sums) < 1e-10),
              info = "Each sample's CLR values must sum to 0")
})

# Check that the correlation matrix is square, symmetric, and has ones on the diagonal
test_that("PCLR is a square symmetric matrix with ones on the diagonal", {
  otu_clr <- CLR(otu_raw + 0.5)
  # Compute the Pearson correlation matrix on CLR-transformed data
  PCLR    <- cor(otu_clr)
  
  # check that rows and column numbers are equal (square matrix check)
  expect_equal(nrow(PCLR), ncol(PCLR))
  # Matrix dimension must be equal the number of OTUs
  expect_equal(nrow(PCLR), ncol(otu_clr))
  # Correlation matrices must be symmetric: cor(i,j) == cor(j,i)
  expect_true(isSymmetric(PCLR))
})

# Check that prevalence values are valid proportions and carry the correct OTU names
test_that("prevalence values are in [0, 1] and correctly named", {
  # Compute fraction of samples where each OTU has count > 0
  otu_prev <- setNames(colSums(otu_raw > 0) / nrow(otu_raw), colnames(otu_raw))
  
  # All prevalence values must be valid proportions
  expect_true(all(otu_prev >= 0 & otu_prev <= 1))
  # Names of the prevalence vector must match OTU column names
  expect_equal(names(otu_prev), colnames(otu_raw))
})



# Check that the info table contains exactly the upper-triangle pairs with no duplicates or self-pairs
test_that("info table contains only upper-triangle pairs (no duplicates, no self-pairs)", {
  otu_clr  <- CLR(otu_raw + 0.5)
  PCLR     <- cor(otu_clr)
  otu_prev <- setNames(colSums(otu_raw > 0) / nrow(otu_raw), colnames(otu_raw))
  
  # Build correlation table, keeping only upper-triangle entries
  info <- data.frame(PCLR) %>%
    rownames_to_column("OTU_I") %>%
    pivot_longer(!OTU_I, names_to = "OTU_J", values_to = "cor") %>%
    # Retain only pairs where OTU_I > OTU_J (upper triangle, no diagonal)
    filter(OTU_I > OTU_J) %>%
    mutate(prev_I = otu_prev[OTU_I],
           prev_J = otu_prev[OTU_J])
  
  # No row should have the same OTU in both columns
  expect_true(all(info$OTU_I != info$OTU_J),
              info = "Self-correlations must be excluded")
  
  # Build a unique key for each pair and check there are no duplicates
  pair_keys <- paste(info$OTU_I, info$OTU_J)
  expect_equal(length(pair_keys), length(unique(pair_keys)),
               label = "Each OTU pair must appear exactly once")
  
  # Upper triangle of an n x n matrix contains exactly n*(n-1)/2 pairs
  n <- ncol(otu_raw)
  expect_equal(nrow(info), n * (n - 1) / 2)
})

# Check that zero-rate columns are valid percentages and complement the prevalence values
test_that("zero_I and zero_J are percentages complementary to prevalence", {
  otu_prev <- setNames(colSums(otu_raw > 0) / nrow(otu_raw), colnames(otu_raw))
  otu_clr  <- CLR(otu_raw + 0.5)
  PCLR     <- cor(otu_clr)
  
  info <- data.frame(PCLR) %>%
    rownames_to_column("OTU_I") %>%
    pivot_longer(!OTU_I, names_to = "OTU_J", values_to = "cor") %>%
    filter(OTU_I > OTU_J) %>%
    mutate(prev_I = otu_prev[OTU_I], prev_J = otu_prev[OTU_J]) %>%
    # Compute zero-rate as percentage of samples where the OTU is absent
    mutate(zero_I = 100 * round(1 - prev_I, 2),
           zero_J = 100 * round(1 - prev_J, 2))
  
  # Zero-rates must be valid percentages in [0, 100]
  expect_true(all(info$zero_I >= 0 & info$zero_I <= 100))
  expect_true(all(info$zero_J >= 0 & info$zero_J <= 100))
  # zero_I must equal 100 * round(1 - prev_I, 2) up to floating point tolerance
  expect_true(all(abs(info$zero_I - 100 * round(1 - info$prev_I, 2)) < 1e-9))
})


# Check that info_filt only retains pairs satisfying both the prevalence and correlation thresholds
test_that("info_filt respects prevalence <= 0.5 and |cor| >= 0.4 thresholds", {
  otu_prev <- setNames(colSums(otu_raw > 0) / nrow(otu_raw), colnames(otu_raw))
  otu_clr  <- CLR(otu_raw + 0.5)
  PCLR     <- cor(otu_clr)
  
  info <- data.frame(PCLR) %>%
    rownames_to_column("OTU_I") %>%
    pivot_longer(!OTU_I, names_to = "OTU_J", values_to = "cor") %>%
    filter(OTU_I > OTU_J) %>%
    mutate(prev_I = otu_prev[OTU_I], prev_J = otu_prev[OTU_J])
  
  # Keep only low-prevalence pairs (<=50%) with strong correlation (|r|>=0.4)
  info_filt <- info %>%
    filter(prev_I <= .5, prev_J <= .5, abs(cor) >= .4)
  
  # All retained pairs must satisfy the prevalence upper bound for both OTUs
  expect_true(all(info_filt$prev_I   <= .5))
  expect_true(all(info_filt$prev_J   <= .5))
  # All retained pairs must have absolute correlation >= 0.4
  expect_true(all(abs(info_filt$cor) >= .4))
})


# Check that the case_when logic assigns the correct label to every detection pattern
test_that("detection labels cover all four cases and are mutually exclusive", {
  # Build a minimal 4-row matrix covering all zero/non-zero combinations
  otu_test <- matrix(
    c(0, 0,   # both absent
      0, 5,   # only OTU_A absent
      3, 0,   # only OTU_B absent
      4, 6),  # both present
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("OTU_A", "OTU_B"))
  )
  
# Classify each row by its detection pattern
  detection <- as_tibble(otu_test) %>%
    mutate(detection = case_when(
      OTU_A == 0 & OTU_B == 0 ~ "Both are 0",
      OTU_A == 0              ~ "OTU_A is 0",
      OTU_B == 0              ~ "OTU_B is 0",
      OTU_A  > 0 & OTU_B  > 0 ~ "Both are > 0"
    ))
  
# Verify each row received the expected label
  expect_equal(detection$detection[1], "Both are 0")
  expect_equal(detection$detection[2], "OTU_A is 0")
  expect_equal(detection$detection[3], "OTU_B is 0")
  expect_equal(detection$detection[4], "Both are > 0")
  # No sample should be left unclassified
  expect_false(any(is.na(detection$detection)),
               info = "Every sample must receive a detection label")
})

# Check that the detection vector length matches the number of samples in the OTU matrix
test_that("detection vector has the same number of rows as the OTU matrix", {

  # Classify each sample based on the first two OTUs in otu_raw
  detection <- as_tibble(otu_raw) %>%
    mutate(detection = case_when(
      OTU_1 == 0 & OTU_2 == 0 ~ "Both are 0",
      OTU_1 == 0              ~ "OTU_1 is 0",
      OTU_2 == 0              ~ "OTU_2 is 0",
      TRUE                    ~ "Both are > 0"
    )) %>%
    select(detection)
  
  # The detection column must have one row per sample
  expect_equal(nrow(detection), nrow(otu_raw))
})