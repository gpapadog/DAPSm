#' DAPSm with optimal weight
#' Function that chooses the optimal weight and fits DAPSm.
#' 
#' @param dataset
#' Data frame including treatment, outcome, coordinates, observed confounders,
#' and propensity score estimates as 'prop.scores'.
#' @param caliper
#' A caliper for the DAPS Score difference of matched pairs. Defaults to 0.1.
#' Scalar.
#' @param coords.cols
#' If the columns of coordinates are not named 'Longitude' and 'Latitude',
#' coords.columns are the column indices corresponding to longitude and latitude
#' accordingly.
#' @param cov.cols
#' If the weight is set to 'optimal', standardized difference of means will be
#' calculated on the columns whose indices are in cov.cols.
#' @param cutoff
#' The cutoff of standardized difference of means under which the covariates are
#' considered balanced. Specify when weight is set to 'optimal'. Defaults to 0.1.
#' @param trt.col The index of the column in the dataset including the binary
#' treatment. Necessary when the column is not named 'X'.
#' @param w_tol
#' Tolerance on the choice of the optimal weight. Only needed when weight is
#' 'optimal'. Defaults to 0.01.
#' @param distance
#' Function te takes in the distance matrix and returns the standardized distance
#' matrix. Defaults to the function that subtracks the minimum and divides by the
#' range.
#' @param caliper_type
#' Whether we want the caliper to be on DAPS or on the PS. caliper_type must
#' either be 'DAPS', or 'PS'.
#' @param quiet
#' Whether we want to print the performance of weights.
#' @param coord_dist
#' Set to true when we want to use a distance function that calculates the
#' spherical distance of points instead of Euclidean. Defaults to FALSE.
#' @param matching_algorithm
#' Argument with options 'optimal', or 'greedy'. The optimal choice uses the optmatch R
#' package to acquire the matches based on propensity score difference and a caliper on
#' distance. The greedy option matches treated and control units sequentially, starting
#' from the ones with the smallest propensity score difference. Defaults to 'optimal'.
#' @param remove.unmatchables Logical. Argument of the optmatch function. Defaults to
#' FALSE. If set to FALSE, the matching fails unless all treated units are matched. If
#' set to TRUE, matching might return matches only for some of the treated units.
#' 
#' @return List of weight chosen, matched dataset, standardized difference of
#' the columns in cov.cols, indices of matched treated and controls.
#' 
#' @examples
#' data('toyData')
#' toyData$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                            data = toyData)$fitted.values
#' daps_opt <- DAPSopt(toyData, caliper = 0.5, coords.cols = c(4, 5),
#'                     cov.cols = 6:9, trt.col = 1, matching_algorithm = 'greedy')
#' class(daps_opt)
#' names(daps_opt)
DAPSopt <- function(dataset, caliper, coords.cols, cov.cols, cutoff = 0.1,
                    trt.col = NULL, w_tol = 0.01, distance = StandDist,
                    caliper_type = c('DAPS', 'PS'), quiet = FALSE,
                    coord_dist = FALSE, matching_algorithm = c('optimal', 'greedy'),
                    remove.unmatchables = FALSE) {
  
  caliper_type <- match.arg(caliper_type)
  matching_algorithm <- match.arg(matching_algorithm)
  
  dataset <- FormDataset(dataset, trt.col = trt.col, out.col = NULL,
                         ignore.cols = NULL)
  # What we will return
  r <- NULL
  
  interval <- c(0, 1)
  
  while ((interval[2] - interval[1]) > (w_tol / 2)) {
    if (!quiet) {
      print(interval)
    }
    x <- WeightChoice(dataset = dataset, caliper = caliper, coords.col = coords.cols,
                      cov.cols = cov.cols, cutoff = cutoff, interval = interval,
                      distance = distance, caliper_type = caliper_type,
                      coord_dist = coord_dist, matching_algorithm = matching_algorithm,
                      remove.unmatchables = remove.unmatchables)
    interval <- x$new_interval
    if (x$success) {
      r$weight <- x$weight
      r$pairs <- x$pairs
      r$stand_diff <- x$stand_diff
      r$ind_trt <- x$ind_trt
      r$ind_cnt <- x$ind_cnt
    }
  }
  if (is.null(r)) {
    warning('Standardized balance not achieved. Weight set to 1.')
    daps.out <- dist.ps(treated = dataset[dataset$X == 1, ],
                        control = dataset[dataset$X == 0, ], caliper = caliper,
                        weight = 1, coords.columns = coords.cols, distance = distance,
                        caliper_type = caliper_type, coord_dist = coord_dist,
                        matching_algorithm = matching_algorithm, 
                        remove.unmatchables = remove.unmatchables)
    r$weight <- 1
    pairs.out        <- daps.out$match
    names(pairs.out) <- rownames(daps.out)
    pairs.out        <- na.omit(pairs.out)
    pairs.daps  <- dataset[c(as.numeric(names(pairs.out)),
                             as.numeric(pairs.out)), ]
    diff_mat <- as.matrix(pairs.daps[, cov.cols])
    stand_diff <- apply(diff_mat, 2, function(x) {
      (mean(x[pairs.daps$X == 1]) - mean(x[pairs.daps$X == 0])) /
        sd(x[pairs.daps$X == 1])
    })
    r$stand_diff <- stand_diff
    r$pairs <- pairs.daps
    r$ind_trt <- as.numeric(names(pairs.out))
    r$ind_cnt <- as.numeric(pairs.out)
  }
  return(r)
}
