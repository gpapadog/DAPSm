#' Standardizing a matrix to be between 0 and 1.
#'
#' Takes in a matrix and subtract the minimum element and divides by the range
#' to return a matrix with elements between 0 and 1.
#'
#' @param x Matrix to be standardized. There should be at least two different
#' entries in x.
#'
#' @return A matrix of the same dimensions as the one given.
#' @export
#' 
#' @examples
#' set.seed(1)
#' D <- matrix(rexp(800, rate = 0.5), 20, 40)
#' Dnew <- StandDist(D)
#' hist(D)
#' hist(Dnew)
StandDist <- function(x) {
  standx <- (x - min(x)) / (max(x) - min(x))
  return(standx)
}


#' Standardizing a matrix to be between 0 and 1 using the empirical CDF.
#'
#' Takes in a matrix and returns a matrix of the same dimension with entries
#' equal to the corresponding empirical cdf value.
#'
#' @param x Matrix to be standardized.
#'
#' @return A matrix of the same dimensions as the one given.
#' @export
#' 
#' @examples
#' set.seed(1)
#' D <- matrix(rexp(800, rate = 0.5), 20, 40)
#' Dnew <- EmpCDF(D)
#' hist(D)
#' hist(Dnew)
EmpCDF <- function(x) {
  y <- stats::ecdf(x)(x)
  y <- matrix(y, nrow = nrow(x), ncol = ncol(x))
  return(y)
}
