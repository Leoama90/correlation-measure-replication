# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# loads the main script to make its variables available in the test environment
source(here::here("biases_parameters_examples", "biases_parameters_examples.R"))

# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

test_that("The object isn't exactly what you were looking for...",{

# test Mean Pielou
# generating a 3x3 diagonal matrix to test the mean_Pielou function
test_matrix <- diag(3) 
expect_lte(mean_Pielou(test_matrix), 1, label = NULL, expected.label = NULL)
  
# test the dimension of the correlation matrixes
expect_shape(cor_D5, dim = c(5, 5))
expect_shape(cor_D30, dim = c(30, 30))
expect_shape(cor_D100, dim = c(100, 100))

# test dimension of the ToyModels
expect_length(toy_D5_P100, 9)
expect_length(toy_D5_P50, 9)
expect_length(toy_D30_P50, 9)

})

