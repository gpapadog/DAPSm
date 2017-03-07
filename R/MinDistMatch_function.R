#' Matching on distance matrix
#'
#' Takes in a matrix and runs through it matching rows to columns with small entries.
#' The smallest one is matched first, and it continues on.
#'
#' @param M       A non-negative matrix (T x C)
#' @param caliper The caliper controls how far apart with respect to their M entry
#'                matched pairs can be. When set to NULL, all possible matches with
#'                finite M entry are taken. If set to a number, no entries are
#'                allowed to be matched if larger than caliper.
#'
#' @return An ordered matrix of the row and column indices that are matched.
#'
#' @export
#' @examples
#' set.seed(1)
#' D <- matrix(rexp(800, rate = 0.5), 20, 40)
#' MinDistMatch(D, caliper = NULL)
#' MinDistMatch(D, caliper = 0.1)
MinDistMatch <- function(M, caliper = NULL) {
  
  num_trt <- nrow(M)
  num_con <- ncol(M)
  rownames(M) <- 1:num_trt
  colnames(M) <- 1:num_con
  
  mat <- NULL
  min_m <- min(M)
  
  # If caliper is set to NULL, all finite matches are allowed.
  if (is.null(caliper)) {
    caliper <- max(as.numeric(M)[!is.infinite(as.numeric(M))]) + 1
  }
  
  while (min(M) <= caliper & !is.null(dim(M))) {
    wh_row <- which(apply(M, 1, function(x) any(x == min_m)))[1]
    wh_col <- which(M[wh_row, ] == min_m)[1]
    mat <- rbind(mat, c(rownames(M)[wh_row], colnames(M)[wh_col]))
    M <- M[ - wh_row, - wh_col]
    min_m <- min(M)
  }
  
  # In the end we're left with a vector, so we need to check for the last matched pair.
  if (is.null(dim(M))) {
    wh_col <- which(M == min(M))
    if (any(M[wh_col] <= caliper)) {
      if (num_trt < num_con) {
        mat <- rbind(mat, c(setdiff(1:num_trt, mat[, 1]), names(M)[wh_col]))
      } else if (num_trt >= num_con) {
        mat <- rbind(mat, c(names(M)[wh_col], setdiff(1:num_con, mat[, 1])))
      }
    }
  }
  
  if (is.null(mat)) {
    stop('No matches were found.')
  }
  mat <- mat[order(mat[, 1]), ]
  mat <- matrix(as.numeric(mat), ncol = 2, nrow = nrow(mat))
  colnames(mat) <- c('Row Index', 'Column Index')
  return(mat)
}