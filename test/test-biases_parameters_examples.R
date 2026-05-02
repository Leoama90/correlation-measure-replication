# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# test the dimension of the correlation matrixes
expect_shape(cor_D5, dim = c(5, 5))
expect_shape(cor_D30, dim = c(30, 30))
expect_shape(cor_D100, dim = c(100, 100))

