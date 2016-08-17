# Author: Georgia Papadogeorgou
# Date: 12/14/2015
# Description: Plotting functions for matching procedures.

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


#' Plotting covariate balance before and after matching.
#' 
#' Function that plots the standardized difference of means before and after matching.
#' 
#' @param dtaBef   Data frame include the data before matching.
#' @param dtaAfter Data frame including only the matched pairs data. Can be set to
#' NULL if we want to plot the standardized difference of means for only one dataset.
#' @param cols     Column indeces of the variables we want to check.
#' @param trt      Column index for treatment variable.
#' @param cutoff   Cutoff for standardized difference. Defaults to 0.1.
#' 
#' @return Plot of standardized difference of means before and after matching.
#' @export
PlotBalance <- function(dtaBef, dtaAfter, cols, trt, cutoff = 0.1) {
  
  stand.means <- CalculateBalance(dtaBef, dtaAfter, cols, trt)
  
  colors <- c("chocolate3", "darkolivegreen4")
  dev.off()
  par(oma = c(1, 5, 0, 1))
  plot(stand.means[1, length(cols):1], y = 1:length(cols), ylab="",
       xlab = "Standardized difference of means", 
       xlim = range(stand.means, cutoff, - cutoff, na.rm = TRUE), axes=F,
       pch=19, cex=1.5, col = colors[1])
  if (!is.null(dtaAfter)) {
    points(stand.means[2, length(cols):1], y = 1:length(cols), pch=17, cex=1.5,
           col=colors[2])
  }
  abline(h=0:length(cols))
  par(las=1)
  axis(2, labels = colnames(stand.means)[length(cols):1],
       at = 1:length(cols), cex=0.5)
  axis(1)
  abline(v = c(cutoff, - cutoff), lty = "dashed", col = "red")
  abline(v = 0)
  legend("bottomleft", legend = c("Before", "After"), pch = c(19,17),
         col = colors)
  
}