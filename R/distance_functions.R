StandDist <- function(x) {
  # Function that takes in a matrix and returns the equivalent
  # matrix with elements between 0 and 1, by subtracting the
  # minimum and dividing by the range.
  #
  # Args:
  #  x: Matrix of any dimensions. There should be at least two
  #     different elements in x. Numeric.
  #
  # Returns:
  #  Matrix of the same dimensions as x with entries between
  #  0 and 1.
  
  standx <- (x - min(x)) / (max(x) - min(x))
  return(standx)
}

EmpCDF <- function(x) {
  # Function that takes in a matrix and returns the empirical
  # cdf of every entry of the matrix.
  #
  # Args:
  #  x: Matrix of any dimensions. There should be at least two
  #     different elements in x. Numeric.
  #
  # Returns:
  #  Matrix of the same dimensions as x with empirical cdf.
  
  y <- stats::ecdf(x)(x)
  y <- matrix(y, nrow = nrow(x), ncol = ncol(x))
  return(y)
}
