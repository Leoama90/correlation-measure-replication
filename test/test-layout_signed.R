# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
library(igraph)

# -- test for errors --------
test_that("there is an error here", {
  
  # load function from script
  source(here("script", "method_comparison", "layout_signed.R"))
  
  # set seed for reproducibility
  set.seed(42)
  # generate non-graph dummy variable
  dum <- runif(20, min = 0, max = 720)
  expect_error(LAYOUT_SIGNED(dum))
  
  # test NULL input
  expect_error(LAYOUT_SIGNED(NULL))
})

# -- test layout_signed function --------
test_that("there is a problem with the output of the script", {
  
  # load function from script
  source(here("script", "method_comparison", "layout_signed.R"))  
  
  # generate dummy graph
  g <- igraph::graph_from_edgelist(
    matrix(c(1,2, 2,3, 3,4, 4,1), ncol=2, byrow=TRUE),
    directed = FALSE
  )
  
  # assign weights to arcs
  igraph::E(g)$weight <- c(1, -1, 0.5, -0.3)
  
  # apply function to dummy graph
  lg <- LAYOUT_SIGNED(g)  
  # count columns in output matrix
  nc <- ncol(LAYOUT_SIGNED(g))
  # verify the function has a matrix output
  expect_true(is.matrix(lg))
  # check there are 2 columns
  expect_equal(nc, 2)
  # check number of rows of the output matrix equals original graph nodes
  expect_true(igraph::vcount(g) == nrow(lg))
  # check the values in matrix are double
  expect_all_true(is.double(lg))
  # check the values in matrix are not integers
  expect_all_false(is.integer(lg))
})

# -- test layout_signed function: edge cases --------
test_that("there is a problem with edge cases", {
  
  # load function from script
  source(here("script", "method_comparison", "layout_signed.R"))
  
  # generate graph with all negative weights
  g_neg <- igraph::make_ring(4)
  igraph::E(g_neg)$weight <- c(-1, -2, -0.5, -3)
  
  # check function returns a matrix even with no positive edges
  expect_true(is.matrix(LAYOUT_SIGNED(g_neg)))
  # check number of rows equals original graph nodes even with no positive edges
  expect_equal(nrow(LAYOUT_SIGNED(g_neg)), igraph::vcount(g_neg))
  
  # generate graph with all positive weights
  g_pos <- igraph::make_ring(4)
  igraph::E(g_pos)$weight <- c(1, 2, 0.5, 3)
  
  # check function returns a matrix with all positive edges
  expect_true(is.matrix(LAYOUT_SIGNED(g_pos)))
  # check number of rows equals original graph nodes with all positive edges
  expect_equal(nrow(LAYOUT_SIGNED(g_pos)), igraph::vcount(g_pos))
  
  # generate graph where one node has no positive edges (isolated node)
  g_iso <- igraph::graph_from_edgelist(
    matrix(c(1,2, 2,3, 3,4), ncol=2, byrow=TRUE),
    directed = FALSE
  )
  # node 4 has no positive edges, node 1 has no positive edges
  igraph::E(g_iso)$weight <- c(-1, 1, -0.5)
  
  # check isolated nodes are preserved in the output (delete.vertices = FALSE)
  expect_equal(nrow(LAYOUT_SIGNED(g_iso)), igraph::vcount(g_iso))
  
  # generate single node graph
  g_single <- igraph::make_empty_graph(n = 1, directed = FALSE)
  
  # check output matrix has 1 row and 2 columns for single node graph
  expect_equal(nrow(LAYOUT_SIGNED(g_single)), 1)
  expect_equal(ncol(LAYOUT_SIGNED(g_single)), 2)
})