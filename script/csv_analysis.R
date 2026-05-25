library(tidyverse)
library(here)

# -------- custom function for preliminary analysis --------


analysis <- function(x){
  
  cat("column number is:", ncol(x),      "\n\n") # get column number
  cat("column names are:", colnames(x),  "\n\n") # get column names
  
  cat("row number is:"    , nrow(x),     "\n\n") # get row number
  if (nrow(x) < 10) { # print row names only if number of rows is < 10
  cat("row names are:"    , rownames(x), "\n\n") # get row names  
  }
  cat("number of vectors in the file is:",  length(x), "\n") # gets number of vectors
  # head(x)
}

