#' Calculate balance of covariates before and after matching.
#' 
#' Function that calculates the standardized difference of means.
#' 
#' @param dtaBef
#' Data frame including the data before matching.
#' @param dtaAfter
#' Data frame including only the matched pairs data. Can be set to NULL if
#' we want to calculate the standardized difference of means for only one
#' dataset.
#' @param cols
#' Column indeces of the variables we want to check. Should be the same in
#' both dtaBef and dtaAfter
#' @param trt
#' Column index for treatment variable.
#' 
#' @return A 2 by length(cols) matrix with the standardized difference of
#' means before and after matching. (Second row is empty if dtaAfter is
#' not given.)
#' 
#' @export
CalculateBalance <- function(dtaBef, dtaAfter = NULL, cols, trt) {
  
  dtaBef <- as.data.frame(dtaBef)
  if (!is.null(dtaAfter)) {
    dtaAfter <- as.data.frame(dtaAfter)
  }
  
  stand.means <- array(NA, dim=c(2, length(cols)))
  # before/after match, covariates
  rownames(stand.means) <- c("Before matching", "After matching")
  colnames(stand.means) <- names(dtaBef)[cols]
  
  for (cc in 1:length(cols)){
    stand.means[1, cc] <- (mean(dtaBef[dtaBef[, trt] == 1, cols[cc]], na.rm = T) -
                             mean(dtaBef[dtaBef[, trt] == 0, cols[cc]], na.rm=T)) /
      sd(dtaBef[dtaBef[, trt] == 1, cols[cc]], na.rm=T)
  }
  if (!is.null(dtaAfter)) {
    for (cc in 1:length(cols)){
      stand.means[2, cc] <- (mean(dtaAfter[dtaAfter[, trt] == 1, cols[cc]], na.rm = T) -
                               mean(dtaAfter[dtaAfter[, trt] == 0, cols[cc]], na.rm=T)) /
        sd(dtaAfter[dtaAfter[, trt] == 1, cols[cc]], na.rm=T)
    }
  }
  return(stand.means)
}
