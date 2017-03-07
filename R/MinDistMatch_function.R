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
  
  mat <- NULL
  if (is.null(caliper)) {
    caliper <- Inf
  }
  
  num_trt <- nrow(M)
  num_con <- ncol(M)
  rownames(M) <- 1:num_trt
  colnames(M) <- 1:num_con
  
  while (nrow(M) > 0) {
    
    # Get the control that is closest to each treated unit.
    mins <- apply(M, 1, min)
    drop_rows <- which(mins > caliper)
    
    # Drop treated units that will not be matched within the caliper.
    if (length(drop_rows) > 0) {
      print(paste0(length(drop_rows), ' units do not have a match with this caliper.'))
      M <- M[- drop_rows, ]
      mins <- mins[- drop_rows]
    }
    
    if (nrow(M) > 0) {
      # Order the matrix such that the treated with the closest control is at the top.
      M <- M[order(mins), ]
      mins <- sort(mins)
      
      # Which control is the closes to each treated. 
      wh_mins <- sapply(1:nrow(M), function(x) which(M[x, ] == mins[x])[1])
      
      # Has the control been seen before.
      duplic <- duplicated(wh_mins)
      stop_matching <- which(duplic)[1] - 1
      
      if (is.na(stop_matching)) {  # No control is used twice.
        
        mat <- rbind(mat, cbind(rownames(M), colnames(M)[wh_mins]))
        M <- matrix(NA, nrow = 0, ncol = 1)
        
      } else {  # Some control is the closest to two or more treated units.

        mat <- rbind(mat, cbind(rownames(M)[1:stop_matching],
                                colnames(M)[wh_mins[1:stop_matching]]))
        
        rows <- rownames(M)[- c(1:stop_matching)]
        cols <- colnames(M)[- wh_mins[1:stop_matching]]
        M <- M[-c(1:stop_matching), -wh_mins[1:stop_matching]]
        
        if (!is.matrix(M)) {
          M <- matrix(M, nrow = 1, ncol = length(M))
          rownames(M) <- rows
          colnames(M) <- cols
        }
      }
    }
  }
  
  if (is.null(mat)) {
    stop('No matches were found.')
  }
  
  mat <- matrix(as.numeric(mat), ncol = 2, nrow = nrow(mat))
  mat <- mat[order(mat[, 1]), ]
  colnames(mat) <- c('Row Index', 'Column Index')
  
  return(mat)
}
