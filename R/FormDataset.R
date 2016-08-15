#' Forming a data set.
#'
#' Takes in a data frame and changes the names of the treatment and outcome
#' columns, and drops unnecessary ones.
#'
#' @param dataset The data frame we want to reform.
#' @param ignore.cols The indices of the columns we want to drop. Defaults
#' to NULL. If set to NULL, no columns are dropped.
#' @param out.col The index of the outcome column which will be renamed to
#' 'Y'. Defaults to NULL. If set to NULL, no column will be renamed.
#' @param trt.col The index of the treatment column which will be renamed
#' to 'X'. Defaults to NULL. If set to NULL, no column will be renamed.
#'
#' @return A data frame with renamed the columns of treatment and outcome,
#' and no longer includes some of the columns.
#' 
#' @examples
#' D <- cbind(Trt = rbinom(100, 1, 1/2), Out = rnorm(100),
#'            C1 = rnorm(100), C2 = rnorm(100))
#' D_new <- FormDataset(D, ignore.cols = 4, out.col = 2, trt.col = 1)
FormDataset <- function(dataset, ignore.cols = NULL,
                        out.col = NULL, trt.col = NULL) {
  # Function that takes a dataset, and reforms it if any of the
  # ignore.cols, out.col, or trt.col is not NULL. Specificially,
  # it drops columns in ignore.cols, renames the outcome column
  # to 'Y', and the treatment column to 'X'.
  #
  # Args:
  #  dataset:     Data frame that we want to reform.
  #  ignore.cols: Indeces of columns that should be dropped from
  #               the data frame. If not specified, no columns
  #               are dropped.
  #  out.col:     If out.col is not NULL, the column of index
  #               out.col will be renamed to 'Y'.
  #  trt.col:     If trt.col is not NULL, the column of index
  #               trt.col will be renamed to 'X'.
  #
  # Returns:
  #  A data frame of same number of rows, after dropping columns
  #  in ignore.cols, and renaming columns out.col, trt.col.
  
  if (!is.null(out.col)) {
    names(dataset)[out.col] <- 'Y'
  }
  if (!is.null(trt.col)) {
    names(dataset)[trt.col] <- 'X'
  }
  if (!is.null(ignore.cols)) {
    dataset <- dataset[, - ignore.cols]
  }
  return(dataset)
}

