# demo_datasummary.R
#
# Purpose:
#   This script shows how datasum() (defined in datasummary.R) works,
#   by building a small dummy tibble with a known structure (including
#   an all-zero column, to show how datasum() handles it) and running
#   datasum() on it, both as a tibble and as a matrix (to also show the
#   determinant and minimum eigenvalue, which are only computed for
#   square numeric matrices, not for tibbles).
#
# Input:
#   - a tibble created inside this script, named "dum_tibble"
#
# Output:
#   - printed column count, row count, mean, median, standard deviation,
#     skewness, number of zeroes and their percentage of dum_tibble
#     (determinant and minimum eigenvalue only when run on its matrix
#     version, since dum_tibble itself is not recognized as a matrix)
#
# Used scripts:
#   - datasummary.R

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)


# -------- load the datasummary.R script --------

# load datasum(), the only function this demo needs
source(
  list.files(
    path = here(),
    pattern = "^datasummary\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- body of the demo --------

# separation line
cat("\n#-------------------------------------------------------#\n")
# UI part where the script explains how the dummy tibble was created
cat("the script will now create a dummy tibble with the following columns:\n
      - column a: numbers from 1 to 5 \n
      - column b: the values in column a + 1 \n
      - column c: the values in column a squared minus 2 \n
      - column d: all values are zeroes (deliberately, to show how ",
    "datasum() reports a fully sparse column) \n
      - column e: values in column b minus values in column a \n")
# separation line
cat("\n#-------------------------------------------------------#\n")

# creating the dum_tibble
dum_tibble <- tibble(a = 1:5, b = a + 1, c = a^2 - 2, d = 0, e = b - a)

# message to announce that dum_tibble will be showed
cat("\nThis is the new dum_tibble\n")
# shows the created dum_tibble
print(dum_tibble)
# separation line
cat("\n#-------------------------------------------------------#\n")

# UI message explaining the first call, on the tibble itself
cat("\nRunning datasum() on dum_tibble directly (determinant and minimum ",
    "eigenvalue will be NA, since a tibble is not recognized as a matrix):\n")
# runs the datasum function on the dum_tibble
datasum(dum_tibble)


# -------- running datasum() on the matrix version, for the determinant --------

# separation line
cat("\n#-------------------------------------------------------#\n")
# UI message explaining the second call, on the matrix version
cat("\nRunning datasum() on as.matrix(dum_tibble) instead, so the",
    "determinant and minimum eigenvalue are also computed. Note both",
    "will be exactly 0, since the all-zero column d makes the matrix",
    "singular by construction:\n")
# runs the datasum function on the matrix version of dum_tibble
datasum(as.matrix(dum_tibble))