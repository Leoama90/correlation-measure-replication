#' Centered Log-Ratio Transformation
#' 
#' @description It calculates the centered log-ratio transformation of X. 
#' The zeroes are replaced with 0.65 (https://doi.org/10.1016/j.chemolab.2021.104248).
#'
#' @param X numeric matrix with all elements greater than or equal to 0.
#' @param mar Integer giving the dimension where the function will be applied;
#' 1 for rows and 2 for columns (default 1).
#' 
CLR <- function(X, mar = 1){
  
  # validate that X is a matrix or vector
  if(!is.matrix(X)) stop("X must be a matrix or a vector")
  
  # validate that X is numeric and contains no negative values
  if(!is.numeric(X) | any(X < 0)) stop("X must be numeric with all elements greater than or equal to 0")
  
  # validate that mar is either 1 (rows) or 2 (columns)
  if(!(mar%in%c(1, 2))) stop("mar has as possible values only 1 and 2.")
  
  # replace zeros with 0.65 to avoid log(0) = -Inf
  X[X == 0] <- 0.65
  
  # compute the geometric mean (in log scale) along the specified margin
  ref <- apply(X, mar, function(x) mean(log(x)) )
  
  # subtract the geometric mean from log(X) to obtain the CLR-transformed matrix
  return(as.matrix(log(X) - ref))
  
}