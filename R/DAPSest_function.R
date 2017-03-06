#' Fitting DAPSm with fixed or optimal weight
#' 
#' Function that fits DAPS by choosing the 'optimal' weight, or using the weight
#' specified and returns the causal effect estimate under a specification of
#' caliper. If the weight is set to 'optimal', the optimal weight is chosen as
#' the smallest weight that satisfies standardized difference less than a cutoff.
#' For computational reasons, we do not fit DAPS on a lot of weights, but perform
#' an iterative algorithm to identify the optimal w.
#' 
#' @param dataset
#' Data frame including treatment, outcome, coordinates, propensity score
#' estimates (named prop.scores) and observed confounders.
#' @param out.col 
#' If outcome column name is not 'Y', out.col should be the index of the outcome
#' column.
#' @param trt.col
#' If treatment is not named 'X', trt.col should be set to the index of the
#' treatment column.
#' @param caliper
#' A caliper for the DAPS Score difference of matched pairs. Defaults to 0.1.
#' @param weight
#' Scalar between 0 and 1 or should be set to 'optimal'. Describes the percent
#' of weight to be given on PS difference. 1 - weight is given to standardized
#' distance. If set to 'optimal', the 'optimal' weight is chosen. Defaults to
#' 'optimal'.
#' @param coords.columns
#' If the columns of coordinates are not named 'Longitude' and 'Latitude',
#' coords.columns are the column indices corresponding to longitude and latitude
#' accordingly.
#' @param pairsRet
#' Whether we want to return the information on the matched pairs. Logical.
#' Defaults to FALSE.
#' @param cov.cols
#' If the weight is set to 'optimal', standardized difference of means will be
#' calculated on the columns whose indices are in cov.cols. If the weight is set
#' to a numeric value, then cov.cols can be left NULL, or it can be used if we
#' want the function to return the standardized difference of means of the columns
#' with indices in cov.cols.
#' @param cutoff
#' The cutoff of standardized difference of means under which the covariates are
#' considered balanced. Defaults to 0.1.
#' @param w_tol
#' Tolerance on the choice of the optimal weight. Only needed when weight is
#' 'optimal'. Defaults to 0.01.
#' @param coord_dist
#' Set to true when we want to use a distance function that calculates the
#' spherical distance of points instead of Euclidean. Defaults to FALSE.
#' @param distance
#' Function that takes in the distance matrix and returns the standardized
#' distance matrix. Defaults to the funcion that subtracks the minimum and
#' divides by the range.
#' @param caliper_type
#' Whether we want the caliper to be on DAPS or on the PS. caliper_type must be
#' either 'DAPS', or 'PS'.
#' @param quiet
#' Whether we want to print the choice of weight in DAPS optimal. Defauls to TRUE.
#' @param true_value
#' Numeric. If provided, an indicator of whether the CI covers the true value is
#' returned.
#' @param matching_algorithm
#' Argument with options 'optimal', or 'greedy'. The optimal choice uses the optmatch R
#' package to acquire the matches based on propensity score difference and a caliper on
#' distance. The greedy option matches treated and control units sequentially, starting
#' from the ones with the smallest propensity score difference. Defaults to 'optimal'.
#' @param remove.unmatchables Logical. Argument of the optmatch function. Defaults to
#' FALSE. If set to FALSE, the matching fails unless all treated units are matched. If
#' set to TRUE, matching might return matches only for some of the treated units.
#' 
#' @return A list including: the estimate of the causal effect, and potential
#' standardized difference of means, optimal weight chosen, information on matched
#' pairs.
#' 
#' @export
#' @examples
#' data('toyData')
#' toyData$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                            data = toyData)$fitted.values
#' daps1 <- DAPSest(toyData, out.col = 2, trt.col = 1, caliper = 0.5,
#'                  weight = 'optimal', coords.columns = c(4, 5),
#'                  pairsRet = TRUE, cov.cols = 6:9, cutoff = 0.1,
#'                  w_tol = 0.001, coord_dist = TRUE, caliper_type = 'DAPS')
#' names(daps1)
#' 
#' # Trying for a different value of the caliper
#' daps2 <- DAPSest(toyData, out.col = 2, trt.col = 1, caliper = 0.1,
#'                  weight = 'optimal', coords.columns = c(4, 5),
#'                  pairsRet = TRUE, cov.cols = 6:9, cutoff = 0.1,
#'                  w_tol = 0.001, coord_dist = TRUE, caliper_type = 'DAPS')
#' names(daps2)
#' daps1$weight
#' daps2$weight
DAPSest <- function(dataset, out.col = NULL, trt.col = NULL, caliper = 0.1,
                    weight = 'optimal', coords.columns = NULL, pairsRet = FALSE,
                    cov.cols = NULL, cutoff = 0.1, w_tol = 0.01, coord_dist = FALSE,
                    distance = StandDist, caliper_type = c('DAPS', 'PS'),
                    quiet = FALSE, true_value = NULL,
                    matching_algorithm = c('optimal', 'greedy'),
                    remove.unmatchables = FALSE) {
  
  matching_algorithm <- match.arg(matching_algorithm)
  caliper_type <- match.arg(caliper_type)
  r <- NULL
  r$weight <- weight
  
  # Naming outcome and treatment as 'Y', 'X'.
  dataset <- as.data.frame(dataset)
  dataset <- FormDataset(dataset, ignore.cols = NULL,
                         out.col = out.col, trt.col = trt.col)
  # Fitting DAPS.
  
  if (is.numeric(weight)) {
    daps.out <- dist.ps(treated = dataset[dataset$X == 1, ],
                        control = dataset[dataset$X == 0, ], caliper = caliper,
                        weight = weight, coords.columns = coords.columns,
                        distance = distance, caliper_type = caliper_type,
                        coord_dist = coord_dist,
                        matching_algorithm = matching_algorithm, 
                        remove.unmatchables = remove.unmatchables)
    
    # If no matches were acheived, return missing values.
    if (nrow(daps.out) == 0) {
      warning(paste0('No matches were acheived for weight = ', weight))
      r$ind_trt <- NA
      r$ind_cnt <- NA
      r$est <- NA
      r$se <- NA
      r$pairs <- matrix(NA, nrow = 0, ncol = 12)
      return(r)
    }
    
    # Getting the matched pairs.
    pairs.out        <- daps.out$match
    names(pairs.out) <- rownames(daps.out)
    pairs.out        <- na.omit(pairs.out)
    pairs.daps  <- dataset[c(as.numeric(names(pairs.out)),
                             as.numeric(pairs.out)), ]

    if (!is.null(cov.cols)) {
      diff_mat <- as.matrix(pairs.daps[, cov.cols])
      stand_diff <- apply(diff_mat, 2, function(x) {
        (mean(x[pairs.daps$X == 1]) - mean(x[pairs.daps$X == 0])) /
          sd(x[pairs.daps$X == 1])
      })
      r$stand_diff <- stand_diff
    }
    r$ind_trt <- as.numeric(names(pairs.out))
    r$ind_cnt <- as.numeric(pairs.out)
    
    
  } else if (weight == 'optimal') {
    daps.opt <- DAPSopt(dataset, caliper = caliper, coords.cols = coords.columns,
                        cov.cols = cov.cols, cutoff = cutoff, w_tol = w_tol,
                        distance = distance, caliper_type = caliper_type,
                        quiet = quiet, coord_dist = coord_dist,
                        matching_algorithm = matching_algorithm,
                        remove.unmatchables = remove.unmatchables)
    pairs.daps <- daps.opt$pairs
    r$weight <- daps.opt$weight
    r$stand_diff <- daps.opt$stand_diff
    r$ind_trt <- daps.opt$ind_trt
    r$ind_cnt <- daps.opt$ind_cnt
  }
  
  lmod <- lm(Y ~ X, data = pairs.daps)
  r$est <- lmod$coef[2]
  r$se <- summary(lmod)$coef[2, 2]
  
  if (!is.null(true_value)) {
    r$cover <- (abs(true_value) - r$est < qnorm(0.975) * r$se)
  }
  
  if (pairsRet) {
    which_cols <- c(which(names(dataset) %in% c('X', 'Y', 'prop.scores')))
    which_cols <- c(which_cols, coords.columns)
    if (is.numeric(weight)) {
      pairs <- pairs.daps[1:(nrow(pairs.daps) / 2), which_cols]
      pairs <- cbind(pairs, pairs.daps[(nrow(pairs.daps) / 2 + 1):
                                         nrow(pairs.daps), which_cols])
      names(pairs) <- paste0(rep(c('Trt.', 'Con.'),
                                 each = length(which_cols)), names(pairs))
      pairs$IDtrt <- as.numeric(names(pairs.out))
      pairs$IDcon <- as.numeric(pairs.out)
      r$pairs <- as.matrix(pairs)[, c(1, 3:6, 8:12, 2, 7)]
    } else {  # If the weight is optimal.
      nmatch <- nrow(daps.opt$pairs) / 2
      mtrt <- daps.opt$pairs[1:nmatch, which_cols]
      mcnt <- daps.opt$pairs[(nmatch + 1) : (2 * nmatch), which_cols]
      pairs <- cbind(mtrt, mcnt)
      names(pairs) <- paste0(rep(c('Trt.', 'Con.'),
                                 each = length(which_cols)), names(pairs))
      pairs$IDtrt <- r$ind_trt
      pairs$IDcnt <- r$ind_cnt
      r$pairs <- pairs[, c(1, 3:6, 8:12, 2, 7)]
    }
  }
  return(r)
}
