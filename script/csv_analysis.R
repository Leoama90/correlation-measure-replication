library(tidyverse)
library(here)

# -------- custom function for preliminary analysis --------


analysis <- function(x){
  
  # get column number
  cat("column number is:", ncol(x),      "\n\n") 
  # get column names
  cat("column names are:", colnames(x),  "\n\n") 
  
  # get row number
  cat("row number is:"    , nrow(x),     "\n\n") 
  # print row names only if number of rows is < 10
  if (nrow(x) < 10) { 
  # get row names  
  cat("row names are:"    , rownames(x), "\n\n")   
  }
  
  # gets number of vectors
  cat("number of vectors in the file is:",  length(x), "\n") 
  # head(x)
}

