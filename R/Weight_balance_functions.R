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
#' @export
#' @examples
#' data(toyData)
#' toyData$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                            data = toyData)$fitted.values
#' bal <- CalcDAPSWeightBalance(toyData, weights = seq(0, 1, length.out = 30),
#'                              cov.cols = 6:9, trt.col = 1,
#'                              coords.columns = c(4, 5), caliper = 0.1)
CalcDAPSWeightBalance <- function(dataset, weights, cov.cols, trt.col = NULL,
                                  out.col = NULL, coords.columns, caliper,
                                  caliper_type = c('DAPS', 'PS'), coord_dist = FALSE,
                                  distance = StandDist,
                                  matching_algorithm = c('optimal', 'greedy')) {
  
  caliper_type <- match.arg(caliper_type)
  matching_algorithm <- match.arg(matching_algorithm)

  if (is.null(trt.col)) {
    trt.col <- which(names(dataset) == 'X')
  }
  if (is.null(out.col)) {
    out.col <- which(names(dataset) == 'Y')
  }
  
  balance <- array(NA, dim = c(length(weights), 2, length(cov.cols)),
                   dimnames = list(weight = round(weights, 2), NULL,
                                   names(dataset)[cov.cols]))
  distance_DAPS <- rep(NA, length(weights))
  num_match_DAPS <- rep(NA, length(weights))
  
  pairs <- NULL
  full_pairs <- NULL
  for (ii in 1:length(weights)) {
    if (ii %% 5 == 0) {
      print(ii)
    }
    A <- DAPSest(dataset, out.col = out.col, trt.col = trt.col,
                 coords.columns = coords.columns, weight = weights[ii],
                 caliper = caliper, pairsRet = TRUE, caliper_type = caliper_type,
                 coord_dist = coord_dist, distance = distance,
                 matching_algorithm = matching_algorithm)
    pairs[[ii]] <- as.numeric(A$pairs[, 9:10])
    full_pairs[[ii]] <- A$pairs
    distance_DAPS[ii] <- mean(fields::rdist(A$pairs[, c(3, 4)], A$pairs[, c(7, 8)]))
    A <- dataset[pairs[[ii]], ]
    num_match_DAPS[ii] <- length(pairs[[ii]])
    balance[ii, , ] <- CalculateBalance(dtaBef = as.data.frame(dataset),
                                        dtaAfter = A, trt = trt.col,
                                        cols = cov.cols)
    rm(A)
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
#' @param leg_cex The size of the legend. Defaults to 1.
#' @param plot_title Overall plot title. Defaults to ''.
#' @param title_cex Size of the title. Defaults to 1.
#' 
#' @export
#' @examples
#' data(toyData)
#' toyData$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                            data = toyData)$fitted.values
#' bal <- CalcDAPSWeightBalance(toyData, weights = seq(0, 1, length.out = 30),
#'                              cov.cols = 6:9, trt.col = 1,
#'                              coords.columns = c(4, 5), caliper = 0.3)
#' PlotWeightBalance(bal$balance, weights = seq(0, 1, length.out = 30),
#'                   cutoff = 0.15)
PlotWeightBalance <- function(balance, full_data = - 3, weights, cutoff,
                              axis_cex = 1, mar = c(4, 4, 2, 8), inset = -0.1,
                              ylimit = NULL, leg_cex = 1, plot_title = '',
                              title_cex = 1) {

  if (is.null(ylimit)) {
    ylimit <- range(c(balance[, 2, ], 0, cutoff, balance[1, 1, ]))
  }
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
         legend = dimnames(balance)[[3]], cex = leg_cex, inset = c(inset, 0))
  title(main = plot_title, cex.main = title_cex)
}



#' Choosing the optimal weight and fitting the corresponding DAPSm.
#'
#' After using CalcDAPSWeightBalance() to calculate the balance of covariates for
#' varying values of w, we can choose the w that acheives the optimal crieterion.
#'
#' @param dataset
#' The dataset that was supplied to CalcDAPSWeightBalance() for calculating
#' balance.
#' @param out.col
#' The index of the outcome column if it is not named 'Y' in the dataset.
#' @param trt.col
#' The index of the treatment column if it is not named 'X'.
#' @param balance
#' A 3-dimensional array including the SDM. First dimension is equal to length
#' of weights, second dimension is equal to two corresponding to before and
#' after matching, and third dimension is the covariates. Returned as an element
#' of the list from the function CalcDAPSWeightBalance().
#' @param cutoff
#' The cutoff that is used for ASDM.
#' @param pairs
#' A list where each element corresponds to a weight. Each element is a vector
#' including the row indices of the dataset that are included in the matched
#' dataset for each weight w. 2nd element of the list returned by
#' CalcDAPSWeightBalance().
#' @param full_pairs
#' A list where each element corresponds to a weight. Includes the basic info
#' about the matched pairs. Returned by CalcDAPSWeightBalance() as full_pairs.
#' Can be left NULL.
#' @param distance_DAPS
#' Numeric of length equal to the number of weights. Mean distance of DAPS
#' matches. Can be left NULL. Or use distance_DAPS of CalcDAPSWeightBalance().
#' @param weights The weights that we used to fit DAPSm.
#' @param true_value A value that we wish to check if it is in the confidence interval.
#' 
#' @return A list of: CE estimate and standard error from a linear model
#' including only the matched pairs for the optimal w, the number of matches,
#' mean distance of pairs if distance_DAPS is specified, balance of observed
#' covariates, the chosen weight, and info on the matched pairs if full_pairs
#' is specified.
#' 
#' @export
#' @examples
#' data(toyData)
#' toyData$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                            data = toyData)$fitted.values
#' bal <- CalcDAPSWeightBalance(toyData, weights = seq(0, 1, length.out = 30),
#'                              cov.cols = 6:9, trt.col = 1,
#'                              coords.columns = c(4, 5), caliper = 0.3)
#' PlotWeightBalance(bal$balance, weights = seq(0, 1, length.out = 30),
#'                   cutoff = 0.15)
#' DAPS <- DAPSchoiceModel(toyData, trt.col = 1, balance = bal$balance,
#'                         cutoff = 0.15, pairs = bal$pairs,
#'                         weights = seq(0, 1, length.out = 30))
#' names(DAPS)
#' DAPS$est
DAPSchoiceModel <- function(dataset, out.col = NULL, trt.col = NULL, balance,
                            cutoff = 0.1, pairs, full_pairs = NULL,
                            distance_DAPS = NULL, weights, true_value = NULL) {
  
  if (is.null(out.col)) {
    out.col <- which(names(dataset) == 'Y')
  }
  if (is.null(trt.col)) {
    trt.col <- which(names(dataset) == 'X')
  }
  out_name <- names(dataset)[out.col]
  trt_name <- names(dataset)[trt.col]
  
  r <- NULL
  bal_ach <- which(apply(balance[, 2, ], 1, function(x) sum(abs(x) > cutoff)) == 0)
  if (length(bal_ach) == 0) {  # If balance was not acheived for any.
    warning('Balance not acheived for any weight.')
    r$est <- NA
    r$cover <- NA
    r$weight <- NA
    r$pairs <- NA
    return(r)
  }
  
  # If balance has been acheived, choose the minimum w.
  wh <- min(bal_ach)
  lmod <- lm(as.formula(paste(out_name, '~', trt_name)), data = dataset[pairs[[wh]], ])
  r$est <- lmod$coef[2]
  r$se <- summary(lmod)$coef[2, 2]
  r$num_match <- length(pairs[[wh]]) / 2
  
  if (!is.null(distance_DAPS)) {
    r$distance <- distance_DAPS[wh]
  }
  r$balance <- balance[wh, , ]
  r$weight <- weights[wh]
  if (!is.null(full_pairs)) {
    r$pairs <- full_pairs[[wh]]
  }
  if (!is.null(true_value)) {
    r$cover <- (abs(true_value - r$est) < qnorm(0.975) * r$se)
  }
  
  return(r)
}



#' Plot the effect estimate as a function of w.
#'
#' Plotting the effect estimate from various fit of DAPSm for varying w. The
#' chosen w will be the only red dot. A loess curve is fit to the effect
#' estimates.
#'
#' @param dataset
#' The dataset that was supplied to CalcDAPSWeightBalance() for calculating
#' balance.
#' @param out.col
#' The index of the outcome column if it is not named 'Y' in the dataset.
#' @param trt.col
#' The index of the treatment column if it is not named 'X'.
#' @param weights
#' The weights that we used to fit DAPSm.
#' @param pairs
#' A list where each element corresponds to a weight. Each element is a vector
#' including the row indices of the dataset that are included in the matched
#' dataset for each weight w. 2nd element of the list returned by
#' CalcDAPSWeightBalance().
#' @param chosen_w
#' The weight value that was chosen by DAPSchoiceModel().
#' 
#' @export
#' @examples
#' data(toyData)
#' toyData$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                            data = toyData)$fitted.values
#' bal <- CalcDAPSWeightBalance(toyData, weights = seq(0, 1, length.out = 30),
#'                              cov.cols = 6:9, trt.col = 1,
#'                              coords.columns = c(4, 5), caliper = 0.3)
#' DAPS <- DAPSchoiceModel(toyData, trt.col = 1, balance = bal$balance,
#'                         cutoff = 0.15, pairs = bal$pairs,
#'                         weights = seq(0, 1, length.out = 30))
#' CE <- DAPSWeightCE(dataset = toyData, trt.col = 1,
#'                    weights = seq(0, 1, length.out = 30), pairs = x$pairs,
#'                    chosen_w = DAPS$weight)
#' CE$plot
DAPSWeightCE <- function(dataset, out.col = NULL, trt.col = NULL, weights,
                         pairs, chosen_w) {
  
  out_name <- 'Y'
  trt_name <- 'X'
  if (!is.null(out.col)) {
    out_name <- names(dataset)[out.col]
  }
  if (!is.null(trt.col)) {
    trt_name <- names(dataset)[trt.col]
  }

  CEweight <- matrix(NA, nrow = length(weights), ncol = 3)
  for (ww in 1:length(weights)) {
    lmod <- lm(as.formula(paste(out_name, '~', trt_name)),
               data = dataset[pairs[[ww]], ])
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
