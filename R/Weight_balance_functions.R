#' DAPSm with extensive search for optimal w.
#'
#' Calculates the balance of covariates as a function of weight when fitting
#' DAPS for multiple values of w.
#'
#' @param dataset
#' Data frame including treatment, outcome, coordinates, propensity score
#' estimates (named prop.scores) and observed confounders.
#' @param weights
#' The weights on which we want to fit DAPS. Vector.
#' @param cov.cols
#' The indices of the columns we want to balance.
#' @param trt.col
#' The index of the binary treatment column. If treatment is named 'X' this can
#' be NULL.
#' @param out.col
#' Can be NULL if the outcome column is named 'Y'. Otherwise, it should the
#' index of the outcome column.
#' @param coords.columns
#' If the columns of coordinates are not named 'Longitude' and 'Latitude',
#' coords.columns are the column indices corresponding to longitude and latitude
#' accordingly.
#' @param caliper
#' The value of the caliper that will be used.
#' @param caliper_type 
#' Whether we want the caliper to be on DAPS or on the PS. caliper_type must be
#' either 'DAPS', or 'PS'.
#' @param coord_dist
#' Set to true when we want to use a distance function that calculates the
#' spherical distance of points instead of Euclidean. Defaults to FALSE.
#' @param distance
#' Function that takes in the distance matrix and returns the standardized
#' distance matrix. Defaults to the funcion that subtracks the minimum and
#' divides by the range.
#' 
#' @return A list including: a 3-dimensional array. Dimensions correspond to
#' weights, before/after matching and covariates. Balance can be plotted using
#'  PlotWeightBalance function. A list of the pairs for the different weights.
#'  
#' @example
#' @export
CalcDAPSWeightBalance <- function(dataset, weights, cov.cols, trt.col = NULL,
                                  out.col = NULL, coords.columns, caliper,
                                  caliper_type = c('DAPS', 'PS'),
                                  coord_dist = FALSE, distance = StandDist) {
  
  caliper_type <- match.arg(caliper_type)
  
  balance <- array(NA, dim = c(length(weights), 2, length(cols.balance)),
                   dimnames = list(weight = round(weights, 2), NULL,
                                   names(dataset)[cols.balance]))
  distance_DAPS <- rep(NA, length(weights))
  num_match_DAPS <- rep(NA, length(weights))
  
  pairs <- NULL
  full_pairs <- NULL
  for (ii in 1:length(weights)) {
    if (ii %% 5 == 0) {
      print(ii)
    }
    A <- DAPSest(dataset, out.col = out.col, trt.col = trt.col,
                 coords.columns = coord.cols,
                 weight = weights[ii], caliper = caliper,
                 pairsRet = TRUE, caliper_type = 'DAPS',
                 coord_dist = coord_dist, distance = distance)
    pairs[[ii]] <- as.numeric(A$pairs[, 9:10])
    full_pairs[[ii]] <- A$pairs
    distance_DAPS[ii] <- mean(rdist(A$pairs[, c(3, 4)], A$pairs[, c(7, 8)]))
    A <- dataset[pairs[[ii]], ]
    num_match_DAPS[ii] <- length(pairs[[ii]])
    balance[ii, , ] <- CalculateBalance(dtaBef = as.data.frame(dataset),
                                        dtaAfter = A, trt = trt.col,
                                        cols = cols.balance)
  }
  return(list(balance = balance, pairs = pairs, distance_DAPS = distance_DAPS,
              num_match_DAPS = num_match_DAPS, full_pairs = full_pairs))
}




#' Plotting balance.
#'
#' Plots balance of the covariates as a function of w and before matching.
#'
#' @param balance
#' A 3-dimensional array including the SDM. First dimension is equal to length
#' of weights, second dimension is equal to two corresponding to before and
#' after matching, and third dimension is the covariates. Returned as an element
#' of the list from the function CalcDAPSWeightBalance().
#' @param full_data
#' The value of the x axis where the full data balance will be plotted. Defaults
#' to - 3.
#' @param weights
#' The vector of weights. Will be used to make the xlab.
#' @param cutoff
#' Vertical lines of cutoff used.
#' @param axis_cex
#' The size of the xaxis. Defaults to 1.
#' @param mar
#' Plot margins. Defaults to c(4, 4, 2, 8).
#' @param inset
#' Inset of the legend Defaults to - 0.1.
#' @param ylimit
#' The limit of the y axis.
#' 
#' @example
#' @export
PlotWeightBalance <- function(balance, full_data = - 3, weights, cutoff,
                              axis_cex = 1, mar = c(4, 4, 2, 8), inset = -0.1,
                              ylimit = NULL) {

  if (is.null(ylimit)) {
    ylimit <- range(c(balance[, 2, ], 0, cutoff, balance[1, 1, ]))
  }
  dev.off()
  par(mar = mar)
  plot(1, type = 'n', xlim = c(full_data, length(weights)), axes = FALSE,
       ylim = ylimit, ylab = 'ASDM', xlab = 'Weight')
  axis(2)
  axis(1, labels = c('Full-Data', round(weights, 2)),
       at = c(full_data, 1:length(weights)), cex.axis = axis_cex)
  
  for (cc in 1:dim(balance)[3]) {
    plot_cov <- c(balance[1, 1, cc], balance[, 2, cc])
    lines(1:length(weights), plot_cov[- 1], col = cc, lwd = 1.5,
          lty = ifelse(cc / 9 >= 1, ifelse(cc / 17 >= 1, 5, 3), 1))
    lines(c(full_data, 1), plot_cov[1:2], col = cc, lty = 2)
    points(full_data, plot_cov[1], pch = 16, col = cc)
  }
  abline(h = c(cutoff, - cutoff), lty = 2, lwd = 3)
  par(xpd = TRUE)
  legend('topright', col = 1:dim(balance)[3],
         lty = c(rep(1, 8), rep(3, 8), rep(5, 8)), lwd = 1.5,
         legend = dimnames(balance)[[3]], cex = 0.7, inset = c(inset, 0))
}




DAPSchoiceModel <- function(balance, cutoff, dataset, pairs, full_pairs,
                            distance_DAPS, num_match_DAPS, out.col) {
  
  out_name <- names(dataset)[out.col]
  
  r <- NULL
  wh <- min(which(apply(balance[, 2, ], 1,
                        function(x) sum(abs(x) > cutoff)) == 0))
  lmod <- lm(as.formula(paste(out_name, '~ SnCR')), data = dataset[pairs[[wh]], ])
  r$est <- lmod$coef[2]
  r$se <- summary(lmod)$coef[2, 2]
  r$distance <- distance_DAPS[wh]
  r$num_match <- num_match_DAPS[wh] / 2
  r$balance <- balance[wh, , ]
  r$weight <- weights[wh]
  r$pairs <- full_pairs[[wh]]
  
  r$balance_plots <- GetBalancePlots(dataset, pairs[[wh]])
  
  return(r)
}



DAPSWeightCE <- function(dataset, weights, pairs, chosen_w, out.col) {
  
  out_name <- names(dataset)[out.col]
  
  CEweight <- matrix(NA, nrow = length(weights), ncol = 3)
  for (ww in 1:length(weights)) {
    lmod <- lm(as.formula(paste(out_name, '~ SnCR')), data = dataset[pairs[[ww]], ])
    CEweight[ww, ] <- summary(lmod)$coef[2, 1] +
      summary(lmod)$coef[2, 2] * 1.96 * c(-1, 0, 1)
  }
  CEweight <- as.data.frame(CEweight)
  CEweight$Weight <- weights
  CEweight$num_match <- sapply(1:length(pairs), function(x) length(pairs[[x]]) / 2)
  
  l <- loess(V2 ~ Weight, data = CEweight)
  
  wh <- which(weights == chosen_w)
  names(CEweight)[1:3] <- c('LB', 'Estimate', 'UB')
  CEweight$col <- 'red'
  CEweight$col[wh] <- 'blue'
  
  g <- ggplot(CEweight, aes(x = Weight, y = Estimate, group = 1, color = col)) +
    geom_pointrange(aes(ymin = LB, ymax = UB),
                    data = CEweight, colour = "black") +
    geom_point(shape = 16, size = 3) +
    ggtitle('Causal Effect estimates with 95% confidence intervals') +
    theme(plot.title = element_text(size = rel(1.4)),
          axis.title = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.2))) +
    ylab('Estimate') + theme(legend.position="none") +
    geom_line(aes(x = weights, y = l$fitted), col = 'black')
  
  return(list(CE = CEweight[, 1:5], plot = g))
}