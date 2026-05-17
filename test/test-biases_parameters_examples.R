# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# loads the main script to make its variables available in the test environment
source(here::here("biases_parameters_examples", "biases_parameters_examples.R"))

# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

test_that("Houston we have a problem!",{

# -------- test Mean Pielou --------

# generates a 3x3 matrix with random numbers to test the function
set.seed(42)
test_matrix <- matrix(runif(9), nrow = 3, ncol = 3) 

# -------- testing that the function gives a result lower than 1 and greater than 0 --------

expect_lt(mean_Pielou(test_matrix), 1, label = NULL, expected.label = NULL) 
expect_gt(mean_Pielou(test_matrix), 0, label = NULL, expected.label = NULL)

# -------- testing that for a matrix full of 1, the mean_Pielou gives 1 as result --------

unitary_matrix <- matrix(rep(1), ncol = 3, nrow = 3)
expect_equal(mean_Pielou(unitary_matrix), 1)

  
# -- test the dimension of the correlation matrixes --------

expect_shape(cor_D5, dim = c(5, 5))
expect_shape(cor_D30, dim = c(30, 30))
expect_shape(cor_D100, dim = c(100, 100))

# -- test dimension of the ToyModels --------

expect_length(toy_D5_P100, 9)
expect_length(toy_D5_P50, 9)
expect_length(toy_D30_P50, 9)

})

