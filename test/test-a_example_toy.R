# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# igraph: used to create and analyze network graphs
# https://cran.r-project.org/web/packages/igraph/index.html
library(igraph)
# extrafont: used to use fonts other than the standard PostScript fonts
# https://cran.r-project.org/web/packages/extrafont/index.html
library(extrafont)
# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# load the function to test
source(here("script", "sparsity_effects", "a_example_toy.R"))

# check that adj is a matrix
test_that("This is not a matrix", {
  expect_true(is.matrix(adj))
})

# check that the adjacency matrix is squared
test_that("The adj matrix is not squared", {
  expect_true(ncol(adj) == nrow(adj))
})

# check that the graph has exactly 40 nodes
test_that("The graph does not have 40 nodes", {
  expect_equal(vcount(g), 40)
})

# check that the graph has exactly 1 edge
test_that("The graph does not have exactly 1 edge", {
  expect_equal(ecount(g), 1)
})

# check that the graph is undirected
test_that("The graph is not undirected", {
  expect_false(is_directed(g))
})

# check that the single edge connects node 5 and node 25
test_that("The edge is not between node 5 and node 25", {
  edge_list <- as_edgelist(g)
  expect_true(any((edge_list[,1] == 5 & edge_list[,2] == 25) |
                    (edge_list[,1] == 25 & edge_list[,2] == 5)))
})

# check that the adjacency matrix has exactly 40 rows and 40 columns
test_that("The adjacency matrix is not 40x40", {
  expect_equal(nrow(adj), 40)
  expect_equal(ncol(adj), 40)
})

# check that the sum of all elements is 2
# (one undirected edge means two symmetric non-zero entries)
test_that("The sum of the adjacency matrix is not 2", {
  expect_equal(sum(adj), 2)
})

# check that all elements outside [5,25] and [25,5] are zero
test_that("The adjacency matrix has non-zero elements outside [5,25] and [25,5]", {
  adj_copy <- adj
  adj_copy[5, 25] <- 0
  adj_copy[25, 5] <- 0
  expect_equal(sum(adj_copy), 0)
})

# check that every node has a size attribute assigned
test_that("Not all nodes have a size attribute", {
  expect_equal(length(V(g)$size), vcount(g))
})

# check that all node sizes are within the expected range [3, 12.5]
test_that("Node sizes are not in the range [3, 12.5]", {
  expect_true(all(V(g)$size >= 3 & V(g)$size <= 12.5))
})

# check that node sizes are reproducible when using the same seed
test_that("Node sizes are not reproducible with the same seed", {
  set.seed(4)
  expected_sizes <- runif(n = 40, min = 3, max = 12.5)
  expect_equal(V(g)$size, expected_sizes)
})

# check that the edge color is set to darkgray
test_that("Edge color is not darkgray", {
  expect_equal(E(g)$color, "darkgray")
})

# check that the edge width is set to 2
test_that("Edge width is not 2", {
  expect_equal(E(g)$width, 2)
})

# check that the edge line type is set to 2 (dashed)
test_that("Edge lty is not 2", {
  expect_equal(E(g)$lty, 2)
})

# check that the output PNG file was created
test_that("The output file was not created", {
  expect_true(file.exists(here("script", "sparsity_effects", "example_netw.png")))
})

# check that the output PNG file is not empty
test_that("The output file is empty", {
  expect_gt(file.info(here("script", "sparsity_effects", "example_netw.png"))$size, 0)
})