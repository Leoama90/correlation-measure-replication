#' Get upper triangular values
#'
#' @description Returns as vector the upper triangular values of a square symmetric
#' matrix
#'
#' @param X numeric symmetric matrix
#'
TRIU <- function(X){
  
  # check X is a matrix
  if(!is.matrix(X)) stop("X must be a matrix or a vector")
  # check X is numeric
  if(!is.numeric(X)) stop("X must be numeric")
  # check X is symmetric
  if(!isSymmetric(X)) stop("X must be symmetric")
  
  # extract upper triangular values (excluding diagonal) as vector
  return(X[upper.tri(X, diag = FALSE)])
}