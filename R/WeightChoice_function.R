#' The central DAPS optimal function
#' 
#' Performs DAPS and checks whether balance of the covariates has been achieved,
#' and chooses the next weight that should be checked.
#' 
#' @param dataset
#' Data frame including treatment, outcome, coordinates, propensity score
#' estimates (named prop.scores) and observed confounders.
#' @param trt.col
#' If the treatment column is not named 'X', set trt.col to the index of the
#' column corresponding to the binary treatment.
#' @param caliper
#' A caliper for the DAPS Score difference of matched pairs. Defaults to 0.1.
#' @param coords.cols
#' If the columns of coordinates are not named 'Longitude' and 'Latitude',
#' coords.columns are the column indices corresponding to longitude and latitude
#' accordingly.
#' @param cov.cols
#' If the weight is set to 'optimal', standardized difference of means will be
#' calculated on the columns whose indices are in cov.cols.
#' @param cutoff
#' The cutoff of standardized difference of means under which the covariates are
#' considered balanced. Defaults to 0.1.
#' @param interval
#' The interval in which we are testing the weight. DAPS is fit in the middle of
#' the interval and depending on whether balance is achieved in the middle, the
#' function chooses the left or right half as the next interval in the iterative
#' procedure.
#' @param distance
#' Function that takes in the distance matrix and returns the standardized
#' distance matrix. Defaults to the funcion that subtracks the minimum and
#' divides by the range.
#' @param caliper_type
#' Whether we want the caliper to be on DAPS or on the PS. caliper_type must be
#' either 'DAPS', or 'PS'.
#' @param coord_dist
#' Set to true when we want to use a distance function that calculates the
#' spherical distance of points instead of Euclidean. Defaults to FALSE.
#' 
#' @return List of next interval, matched dataset, standardized difference of
#' the columns in cov.cols, indices of matched treated and controls, whether
#' balance was achieved, and the next interval in the iterative algorithm.
#' 
#' @examples 
#' data('toyData')
#' toyData$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                            data = toyData)$fitted.values
#' r <- WeightChoice(toyData, trt.col = 1, caliper = 0.5, coords.cols = c(4, 5),
#'                   cov.cols= 6:9, cutoff = 0.1, interval = c(0.5, 1),
#'                   distance = StandDist, caliper_type = 'DAPS',
#'                   coord_dist = FALSE)
#' names(r)
WeightChoice <- function(dataset, trt.col = NULL, caliper, coords.cols, cov.cols,
                         cutoff, interval, distance = StandDist,
                         caliper_type, coord_dist = FALSE) {
  
  # I dont need to check caliper_type for this function, since it will not be
  # exported, and it's only used within DAPSopt.
  
  dataset <- FormDataset(dataset, trt.col = trt.col)
  
  r <- NULL
  
  stand_diff <- rep(cutoff + 1, length(cov.cols))
  names(stand_diff) <- colnames(dataset)[cov.cols]
  
  weight <- mean(interval)
  r$weight <- weight
  
  daps.out <- dist.ps(treated = dataset[dataset$X == 1, ],
                      control = dataset[dataset$X == 0, ],
                      caliper = caliper, weight = weight,
                      coords.columns = coords.cols,
                      distance = distance,
                      caliper_type = caliper_type,
                      coord_dist = coord_dist)
  pairs.out        <- daps.out$match
  names(pairs.out) <- rownames(daps.out)
  pairs.out        <- na.omit(pairs.out)
  pairs.daps  <- dataset[c(as.numeric(names(pairs.out)),
                           as.numeric(pairs.out)), ]
  trt <- which(pairs.daps$X == 1)
  cnt <- which(pairs.daps$X == 0)
  if (length(cov.cols) == 1) {
    diff_mat <- pairs.daps[, cov.cols]
    stand_diff <- mean(diff_mat[trt]) - mean(diff_mat[cnt])
    stand_diff <- stand_diff / sd(diff_mat[trt])
  } else {
    stand_diff <- apply(pairs.daps[trt, cov.cols], 2, mean) -
      apply(pairs.daps[cnt, cov.cols], 2, mean)
    stand_diff <- stand_diff / apply(pairs.daps[trt, cov.cols], 2, sd)
  }
  r$stand_diff <- stand_diff
  r$ind_trt <- as.numeric(names(pairs.out))
  r$ind_cnt <- as.numeric(pairs.out)
  r$pairs <- pairs.daps
  r$success <- FALSE
  
  
  if (!any(abs(stand_diff) > cutoff)) {  # If none is above the cutoff.
    r$success <- TRUE
    r$new_interval <- c(interval[1], weight)
  } else {
    r$success <- FALSE
    r$new_interval <- c(weight, interval[2])
  }
  return(r)
}