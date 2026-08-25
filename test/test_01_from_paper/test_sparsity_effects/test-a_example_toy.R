# test-a_example_toy.R
#
# Purpose:
#   this script tests the main features of the a_example_toy.R script,
#   checking:
#       - the adjacency matrix has the expected dimensions and structure
#       - the graph has the expected number of nodes and edges and is undirected
#       - the graph edge connects the expected nodes and the adjacency matrix
#         contains no other edges
#       - node and edge attributes are correctly assigned and reproducible
#       - the expected PNG output file is created and is not empty
#
# Inputs:
#   - a_example_toy.R script (sourced below)
#   
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- load the function to test --------

source(
  list.files(
    path = here(),
    pattern = "^a_example_toy\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- testing adjancency matrix --------

# check that the adjacency matrix has exactly 40 rows and 40 columns
test_that("The adjacency matrix has dimensions 40 x 40", {
  expect_equal(nrow(adj), 40)
  expect_equal(ncol(adj), 40)
})

# check that adj is a matrix
test_that("The adjacency matrix is a matrix", {
  expect_true(is.matrix(adj))
})

# check that the adjacency matrix is squared
test_that("The adjacency matrix is square", {
  expect_true(ncol(adj) == nrow(adj))
})


# -------- test for the graph --------

# check that the graph has exactly 40 nodes
test_that("The graph contains 40 nodes", {
  expect_equal(vcount(g), 40)
})

# check that the graph has exactly 1 edge
test_that("The graph contains exactly 1 edge", {
  expect_equal(ecount(g), 1)
})

# check that the graph is undirected
test_that("The graph is undirected", {
  expect_false(is_directed(g))
})

# check that the single edge connects node 5 and node 25
test_that("The single edge connects nodes 5 and 25", {
  edge_list <- as_edgelist(g)
  expect_true(
    any(
      (edge_list[, 1] == 5 & edge_list[, 2] == 25) |
        (edge_list[, 1] == 25 & edge_list[, 2] == 5)
    )
  )
})

# check that the sum of all elements is 2
# (one undirected edge means two symmetric non-zero entries)
test_that("The adjacency matrix contains exactly one undirected edge", {
  expect_equal(sum(adj), 2)
})

# check that all elements outside [5,25] and [25,5] are zero
test_that("The adjacency matrix contains no other edges", {
  adj_copy <- adj
  adj_copy[5, 25] <- 0
  adj_copy[25, 5] <- 0
  expect_equal(sum(adj_copy), 0)
})

# check that every node has a size attribute assigned
test_that("All graph nodes have a size attribute", {
  expect_equal(length(V(g)$size), vcount(g))
})

# check that all node sizes are within the expected range [3, 12.5]
test_that("All node sizes are within the range [3, 12.5]", {
  expect_true(all(V(g)$size >= 3 & V(g)$size <= 12.5))
})

# check that node sizes are reproducible when using the same seed
test_that("Node sizes are reproducible with the expected seed", {
  set.seed(4)
  expected_sizes <- runif(n = 40, min = 3, max = 12.5)
  expect_equal(V(g)$size, expected_sizes)
})


# -------- test for the edges --------

# check that the edge color is set to darkgray
test_that("The edge color is darkgray", {
  expect_equal(E(g)$color, "darkgray")
})

# check that the edge width is set to 2
test_that("The edge width is 2", {
  expect_equal(E(g)$width, 2)
})

# check that the edge line type is set to 2 (dashed)
test_that("The edge line type is 2", {
  expect_equal(E(g)$lty, 2)
})


# -------- test for the final output of the script --------

# check that the output PNG file was created
test_that("The expected PNG output file exists", {
  expect_true(file.exists(here("Plots", "example_netw.png")))
})

# check that the output PNG file is not empty
test_that("The PNG output file is not empty", {
  expect_gt(file.info(here("Plots", "example_netw.png"))$size, 0)
})