# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# Matrix: needed for nearPD() and eigenvalue checks
# https://cran.r-project.org/web/packages/Matrix/index.html
library(Matrix)


# -- Test Parameters --------

Dimension <- 10  # small dimension for fast testing
n_taxa_connected <- 5
corr_val <- 0.7
total_possible_edges <- Dimension * (Dimension - 1) / 2


# -- Block Correlation Matrix Construction --------

test_that("block correlation matrix is built correctly", {
  
  corM <- diag(1, Dimension, Dimension)
  corM[1:n_taxa_connected, 1:n_taxa_connected] <- corr_val
  diag(corM) <- 1  # restore diagonal after block assignment
  
  # diagonal must be all ones
  expect_true(all(diag(corM) == 1))
  
  # off-diagonal block must match corr_val
  expect_true(all(corM[1:n_taxa_connected, 1:n_taxa_connected][lower.tri(corM[1:n_taxa_connected, 1:n_taxa_connected])] == corr_val))
  
  # outside the block, off-diagonal entries must be zero
  expect_true(all(corM[(n_taxa_connected + 1):Dimension, (n_taxa_connected + 1):Dimension][lower.tri(corM[(n_taxa_connected + 1):Dimension, (n_taxa_connected + 1):Dimension])] == 0))
  
})


# -- Edge Density Calculation --------

test_that("edge density is computed correctly and lies in [0, 1]", {
  
  corM <- diag(1, Dimension, Dimension)
  corM[1:n_taxa_connected, 1:n_taxa_connected] <- corr_val
  diag(corM) <- 1
  
  num_edges <- sum(corM[lower.tri(corM)] != 0)
  edge_density <- num_edges / total_possible_edges
  
  # density must be a single numeric value in [0, 1]
  expect_length(edge_density, 1)
  expect_gte(edge_density, 0)
  expect_lte(edge_density, 1)
  
  # a fully connected matrix must have density 1
  full_corM <- matrix(corr_val, Dimension, Dimension)
  diag(full_corM) <- 1
  full_edges <- sum(full_corM[lower.tri(full_corM)] != 0)
  expect_equal(full_edges / total_possible_edges, 1)
  
})


# -- Positive Semidefiniteness after nearPD --------

test_that("nearPD produces a valid positive semidefinite correlation matrix", {
  
  corM <- diag(1, Dimension, Dimension)
  corM[1:n_taxa_connected, 1:n_taxa_connected] <- corr_val
  diag(corM) <- 1
  
  corM_PD <- nearPD(corM, corr = TRUE, keepDiag = TRUE)$mat
  corM_PD <- as.matrix(corM_PD)
  
  # all eigenvalues must be non-negative
  eigenvalues <- eigen(corM_PD, symmetric = TRUE)$values
  expect_true(all(eigenvalues >= -1e-8))  # small tolerance for floating point
  
  # diagonal must still be all ones
  expect_true(all(abs(diag(corM_PD) - 1) < 1e-8))
  
  # matrix must be symmetric
  expect_true(isSymmetric(corM_PD, tol = 1e-8))
  
})


# -- Gaussian Random Matrix Construction --------

test_that("gaussian random matrix is symmetric and positive semidefinite after correction", {
  
  set.seed(42)
  n <- Dimension
  
  random_matrix <- matrix(rnorm(n * n), n, n)
  symmetric_matrix <- (random_matrix + t(random_matrix)) / 2
  diag(symmetric_matrix) <- 1
  
  # matrix must be symmetric before eigenvalue correction
  expect_true(isSymmetric(symmetric_matrix, tol = 1e-8))
  
  # apply eigenvalue correction
  eigen_values <- eigen(symmetric_matrix)
  eigen_values$values[eigen_values$values < 0] <- 1e-10
  correlation_matrix <- eigen_values$vectors %*% diag(eigen_values$values) %*% t(eigen_values$vectors)
  
  correlation_matrix <- nearPD(correlation_matrix, corr = TRUE, keepDiag = TRUE)$mat
  correlation_matrix <- as.matrix(correlation_matrix)
  
  # all eigenvalues must be non-negative after correction
  eigenvalues <- eigen(correlation_matrix, symmetric = TRUE)$values
  expect_true(all(eigenvalues >= -1e-8))
  
  # diagonal must be all ones
  expect_true(all(abs(diag(correlation_matrix) - 1) < 1e-8))
  
})