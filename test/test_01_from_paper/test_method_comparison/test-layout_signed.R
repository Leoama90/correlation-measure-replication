# test-layout_signed.R
#
# Purpose:
#   this script tests the most important features of the layout_signed.R script,
#   checking:
#       - the function gives an error when receiving invalid inputs
#       - the function returns a valid matrix with the expected dimensions
#       - the output preserves the number of nodes for graphs with positive,
#         negative, or mixed edge weights
#       - the function correctly handles edge cases, including isolated nodes
#         and single-node graphs
#
# Inputs:
#   - layout_signed.R script (sourced below)
#   - dummy graphs generated with igraph for testing
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# igraph: graph analysis for R
# https://cran.r-project.org/web/packages/igraph/index.html
library(igraph)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- load the function to test --------

source(
  list.files(
    path = here(),
    pattern = "^layout_signed\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- test for errors --------

test_that("LAYOUT_SIGNED returns an error for invalid inputs", {
  # set seed for reproducibility
  set.seed(42)

  # generate non-graph dummy variable
  dum <- runif(20, min = 0, max = 720)

  # LAYOUT_SIGNED should reject non-graph input
  expect_error(LAYOUT_SIGNED(dum))

  # test NULL input
  expect_error(LAYOUT_SIGNED(NULL))
})


# -------- test layout_signed function --------

test_that("LAYOUT_SIGNED returns a valid 2D coordinate matrix", {
  # generate dummy graph
  g <- graph_from_edgelist(
    matrix(c(1, 2, 2, 3, 3, 4, 4, 1), ncol = 2, byrow = TRUE),
    directed = FALSE
  )

  # assign weights to arcs
  E(g)$weight <- c(1, -1, 0.5, -0.3)

  # apply function to dummy graph
  lg <- LAYOUT_SIGNED(g)

  # count columns in output matrix
  nc <- ncol(lg)

  # verify the function has a matrix output
  expect_true(is.matrix(lg))

  # check there are 2 columns
  expect_equal(nc, 2)

  # check number of rows of the output matrix equals original graph nodes
  expect_true(vcount(g) == nrow(lg))

  # check the values in matrix are double
  expect_all_true(is.double(lg))

  # check the values in matrix are not integers
  expect_all_false(is.integer(lg))
})


# -------- test layout_signed function: edge cases --------

test_that("LAYOUT_SIGNED handles edge cases correctly", {
  # generate graph with all negative weights
  g_neg <- make_ring(4)
  E(g_neg)$weight <- c(-1, -2, -0.5, -3)

  # check function returns a matrix even with no positive edges
  expect_true(is.matrix(LAYOUT_SIGNED(g_neg)))

  # check number of rows equals original graph nodes even with no positive edges
  expect_equal(nrow(LAYOUT_SIGNED(g_neg)), vcount(g_neg))

  # generate graph with all positive weights
  g_pos <- make_ring(4)
  E(g_pos)$weight <- c(1, 2, 0.5, 3)

  # check function returns a matrix with all positive edges
  expect_true(is.matrix(LAYOUT_SIGNED(g_pos)))

  # check number of rows equals original graph nodes with all positive edges
  expect_equal(nrow(LAYOUT_SIGNED(g_pos)), vcount(g_pos))

  # generate graph where one node has no positive edges (isolated node)
  g_iso <- graph_from_edgelist(
    matrix(c(1, 2, 2, 3, 3, 4), ncol = 2, byrow = TRUE),
    directed = FALSE
  )

  # node 4 has no positive edges, node 1 has no positive edges
  E(g_iso)$weight <- c(-1, 1, -0.5)

  # check isolated nodes are preserved in the output (delete.vertices = FALSE)
  expect_equal(nrow(LAYOUT_SIGNED(g_iso)), vcount(g_iso))

  # generate single node graph
  g_single <- make_empty_graph(n = 1, directed = FALSE)

  # check output matrix has 1 row and 2 columns for single node graph
  expect_equal(nrow(LAYOUT_SIGNED(g_single)), 1)
  expect_equal(ncol(LAYOUT_SIGNED(g_single)), 2)
})
