# demo_pseudocounts.R
#
# Purpose:
#   The aim of this script is to show how pseudocount.R works, printing
#   on screen information about the function's behaviour before running
#   it on a small dummy OTU table, checking:
#       - how raw counts are converted to row-wise proportions
#       - how the row-specific detection limit (1/library size) is used
#       - how zeroes are replaced with a pseudocount, based on the
#         percentage threshold entered by the user at the prompt
#       - that non-zero entries stay unchanged as proportions
#
# Inputs:
#   - pseudocount.R (sourced below)
#   - a threshold value entered interactively at the prompt (pseudocount()
#     asks for it directly; this demo does not set one itself)
#
# Outputs:
#   - a series of information printed on screen, and a demo dataset
#     (demo_result) left in the environment for further inspection

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -------- load the function to demonstrate --------

source(
  list.files(
    path = here(),
    pattern = "^pseudocount\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- descriptive part --------

cat("\n#------------------------------------------------------------------------------#\n")
cat("\npseudocount() converts raw counts to row-wise proportions, then\n")
cat("replaces every zero with a small, row-specific pseudocount, so that\n")
cat("later steps like CLR (which requires strictly positive values) can\n")
cat("be applied without errors.\n")
cat("#------------------------------------------------------------------------------#\n")


# -------- parameters explanation --------

cat("\npseudocount() takes a single argument, x: a numeric matrix/data\n")
cat("frame of non-negative counts (samples on rows, taxa on columns).\n")
cat("It then asks, interactively, for a threshold between 0 and 1: each\n")
cat("zero is replaced with that fraction of the sample's detection limit\n")
cat("(1 / library size). E.g. entering 0.5 replaces a zero with half of\n")
cat("the smallest proportion technically detectable in that sample.\n")


# -------- demo data --------

cat("\n#------------------------------------------------------------------------------#\n")
cat("\nBuilding a small demo OTU table (3 samples x 3 OTUs), with one\n")
cat("zero in each row, in a different column each time:\n\n")

demo_otu <- matrix(
  c(
    10, 0, 5,
     2, 3, 0,
     0, 1, 20  
  ), ncol = 3,
  nrow = 3, byrow = TRUE,
  dimnames = list(NULL, c("OTU1", "OTU2", "OTU3"))
)
print(demo_otu)

cat("\nLibrary sizes (row sums):", rowSums(demo_otu), "\n")
cat("Detection limits (1 / library size):", round(1 / rowSums(demo_otu), 4), "\n")


# -------- run pseudocount() interactively --------

cat("\n#------------------------------------------------------------------------------#\n")
cat("\nCalling pseudocount(demo_otu) now. You will be asked to enter a\n")
cat("threshold between 0 and 1 (try 0.5 for a first look).\n")
cat("#------------------------------------------------------------------------------#\n\n")

demo_result <- pseudocount(demo_otu)


# -------- show the result --------

cat("\n#------------------------------------------------------------------------------#\n")
cat("\nResult (demo_result): raw counts converted to proportions, with\n")
cat("zeroes replaced by their row-specific pseudocount:\n\n")
print(demo_result)

cat("\nNotice that non-zero entries are simple proportions (e.g. 10/15\n")
cat("in row 1), while the former zero in each row now holds a small\n")
cat("positive value instead, computed from the threshold you entered.\n")