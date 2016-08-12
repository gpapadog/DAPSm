#' Linear regression
#'
#' Runs an OLS regression not unlike \code{\link{lm}}
#'
#' @param y response vector (1 x n)
#' @param X covariate matrix (p x n) with no intercept
#'
#' @return A list with 4 elements: coefficients, vcov, sigma, df
#'
#' @examples
#' data(mtcars)
#' X <- as.matrix(mtcars[, c("cyl", "disp", "hp")])
#' y <- mtcars[, "mpg"]
#' linreg(y, X)
#'
#' @export
#'
#'

MinDistMatch <- function(M, caliper = NULL) {
  # Function that takes in a matrix M and runs through it matching observations
  # with small value of M. The smallest one is matched first, and it continues in
  # such way.
  #
  # Args:
  #
  #  M:       A non-negative matrix.
  #  caliper: The caliper controls how far apart with respect to their M
  #           entry matched pairs can be. When set to NULL, all possible
  #           matches with finite M entry are taken.
  #
  # Returns:
  #  An ordered matrix of the row and column indeces that are matched.
  
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
  
  # In the end we are left with a vector, so we need to check for
  # the last matched pair.
  if (is.null(dim(M))) {
    wh_col <- which(M == min(M))
    if (M[wh_col] < caliper) {
      if (num_trt < num_con) {
        mat <- rbind(mat, c(setdiff(1:num_trt, mat[, 1]),
                            names(M)[wh_col]))
      } else if (num_trt >= num_con) {
        mat <- rbind(mat, c(names(M)[wh_col],
                            setdiff(1:num_con, mat[, 1])))
      }
    }
  }
  
  if (is.null(mat)) {
    stop('No matches were found.')
  }
  mat <- mat[order(mat[, 1]), ]
  mat <- matrix(as.numeric(mat), ncol = 2, nrow = nrow(mat))
  return(mat)
}





dist.ps <- function(treated, control, caliper = 0.1, weight = 0.8,
                    coords.columns = NULL, distance = StandDist,
                    caliper_type, coord_dist = FALSE) {
  # Function that takes in a two data frames one with treatment and one
  # with control units including the variables: Longitude, Latitude and
  # propensity scores (as prop.scores) and returns a matrix of the
  # matched pairs using DAPS, and information on PS difference, matching
  # difference, distance. Weight on PS can be specified, as well as DAPS
  # caliper.
  #
  # Args:
  #  treated:        A data frame include the treated units and the
  #                  variables: 'Longitude', 'Latitude' and propensity
  #                  scores (named 'prop.scores'). The rownames of treated
  #                  should be the unit ids.
  #  control:        Control units. Same variables as in treated.
  #  caliper:        A caliper of DAPS for matching.
  #  weight:         Number between 0 and 1, percentage of matching weight
  #                  to be given on propensity score difference.
  #  coords.columns: If the columns of coordinates are not named
  #                  'Longitude', 'Latitude', coords.cols should be the
  #                  column indeces corresponding to longitude and
  #                  latitude accordingly.
  #  distance:       Function that takes in the distance matrix and returns
  #                  the standardized distance matrix. Defaults to the
  #                  funcion that subtracks the minimum and divides by
  #                  the range.
  #  caliper_type:   Whether we want the caliper to be on DAPS or on the PS.
  #                  caliper_type must either be 'DAPS', or 'PS'.
  #  coord_dist:     Set to true when we want to use a distance function that
  #                  calculates the spherical distance of points instead of
  #                  euclidean. Defaults to FALSE.
  #
  # Returns:
  #  A dataframe, where each row corresponds to each treated unit, and
  #  includes: the control unit to which it was matched, their propensity
  #  score difference, their DAPS difference, their distance, their
  #  standardized distance.
  
  require(fields)  # In order to use rdist().
  
  if (!is.null(coords.columns)) {
    names(treated)[coords.columns] <- c('Longitude', 'Latitude')
    names(control)[coords.columns] <- c('Longitude', 'Latitude')
  }
  
  if (coord_dist) {  # Spherical distance on the globe.
    dist.mat <- rdist.earth(cbind(treated$Longitude, treated$Latitude),
                            cbind(control$Longitude, control$Latitude))
  } else {
    dist.mat <- rdist(cbind(treated$Longitude, treated$Latitude),
                      cbind(control$Longitude, control$Latitude))
  }
  stand.dist.mat <- distance(dist.mat)
  ps.diff <- t(sapply(1:nrow(treated),
                      function(x) {
                        abs(treated$prop.scores[x] -
                              control$prop.scores)
                      }))
  dapscore <- (1 - weight) * stand.dist.mat + weight * ps.diff
  
  # Creating the matrix we will use for matching, depending on
  # whether the caliper is set on DAPS or on the PS.
  if (caliper_type == 'DAPS') {
    caliper <- caliper * sd(dapscore)
    M <- ifelse(dapscore <= caliper, dapscore, Inf)
  } else if (caliper_type == 'PS') {
    caliper <- caliper * sd(c(treated$prop.scores, control$prop.scores))
    M <- ifelse(ps.diff <= caliper, ps.diff, Inf)
    M <- (1 - weight) * stand.dist.mat + weight * M
  } else {
    stop('Caliper type should be set to DAPS or PS.')
  }
  
  pairs <- MinDistMatch(M, caliper = NULL)
  matched_trt <- pairs[, 1]
  matched_con <- pairs[, 2]
  
  # Where we will save the results.
  mat <- data.frame(match = rep(NA, dim(treated)[1]),
                    distance = rep(NA, dim(treated)[1]),
                    prop.diff = rep(NA, dim(treated)[1]),
                    match.diff = rep(NA, dim(treated)[1]),
                    stand.distance = rep(NA, dim(treated)[1]))
  rownames(mat) <- rownames(treated)
  
  mat$match[matched_trt] <- rownames(control)[matched_con]
  for (ii in 1:length(matched_trt)) {
    wh_trt <- matched_trt[ii]
    wh_con <- matched_con[ii]
    mat$stand.distance[wh_trt] <- stand.dist.mat[wh_trt, wh_con]
    mat$prop.diff[wh_trt] <- treated$prop.scores[wh_trt] -
      control$prop.scores[wh_con]
    mat$match.diff[wh_trt] <- dapscore[wh_trt, wh_con]
    mat$distance[wh_trt] <- dist.mat[wh_trt, wh_con]
  }
  
  return(mat)
}





# The only thing that changed in DAPSopt and WeightChoice is that we
# call NEW_dist.ps instead of DAPS, and we can drop the argument perm.

DAPSopt <- function(dataset, caliper, coords.cols, cov.cols, cutoff = 0.1,
                    w_tol = 0.05, distance = StandDist, caliper_type,
                    quiet = FALSE, coord_dist = FALSE) {
  # Function that chooses the optimal weight and fits DAPS.
  #
  # Args:
  #  dataset:       Data frame including treatment, outcome, coordinates,
  #                 and observed confounders.
  #  caliper:       A caliper for the DAPS Score difference of matched
  #                 pairs. Defaults to 0.1. Scalar.
  #  coords.cols:   If the columns of coordinates are not named 'Longitude'
  #                 and 'Latitude', coords.columns are the column indeces
  #                 corresponding to longitude and latitude accordingly.
  #  cov.cols:      If the weight is set to 'optimal', standardized
  #                 difference of means will be calculated on the columns
  #                 whose indeces are in cov.cols.
  #  cutoff:        The cutoff of standardized difference of means under
  #                 which the covariates are considered balanced. Specify
  #                 when weight is set to 'optimal'. Defaults to 0.1.
  #  w_tol:         Tolerance on the choice of the optimal weight. Only
  #                 needed when weight is 'optimal'. Defaults to 0.05.
  #  distance:      Function te takes in the distance matrix and returns
  #                 the standardized distance matrix. Defaults to the
  #                 function that subtracks the minimum and divides by
  #                 the range.
  #  caliper_type:  Whether we want the caliper to be on DAPS or on the PS.
  #                 caliper_type must either be 'DAPS', or 'PS'.
  #  quiet:         Whether we want to print the performance of weights.
  #  coord_dist:    Set to true when we want to use a distance function that
  #                 calculates the spherical distance of points instead of
  #                 euclidean. Defaults to FALSE.
  #
  # Returns:
  #  List of weight chosen, matched dataset, standardized difference of
  #  the columns in cov.cols, indeces of matched treated and controls.
  
  
  # What we will return
  r <- NULL
  
  interval <- c(0, 1)
  
  while ((interval[2] - interval[1]) > (w_tol * 2)) {
    if (!quiet) {
      print(interval)
    }
    x <- WeightChoice(dataset, caliper, coords.cols, cov.cols,
                      cutoff, interval, distance = distance,
                      caliper_type, coord_dist = coord_dist)
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
                        control = dataset[dataset$X == 0, ],
                        caliper = caliper, weight = 1,
                        coords.columns = coords.cols,
                        distance = distance,
                        caliper_type = caliper_type,
                        coord_dist = coord_dist)
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


WeightChoice <- function(dataset, caliper, coords.cols, cov.cols,
                         cutoff, interval, distance = StandDist,
                         caliper_type, coord_dist = FALSE) {
  # Function that performs DAPS and checks whether balance of the
  # covariates has been achieved, and chooses the next weight that
  # should be checked.
  #
  # Args:
  #  dataset:      Data frame including treatment, outcome, coordinates,
  #                and observed confounders.
  #  caliper:      A caliper for the DAPS Score difference of matched
  #                pairs. Defaults to 0.1. Scalar.
  #  coords.cols:  If the columns of coordinates are not named 'Longitude'
  #                and 'Latitude', coords.columns are the column indeces
  #                corresponding to longitude and latitude accordingly.
  #  cov.cols:     If the weight is set to 'optimal', standardized
  #                difference of means will be calculated on the columns
  #                whose indeces are in cov.cols.
  #  cutoff:       The cutoff of standardized difference of means under
  #                which the covariates are considered balanced. Specify
  #                when weight is set to 'optimal'. Defaults to 0.1.
  #  interval:     The interval in which we are testing the weight. DAPS
  #                is fit in the middle of the interval and depending on
  #                whether balance is achieved in the middle, the
  #                function chooses the left or right half as the next
  #                interval in the iterative procedure.
  #  distance:     Function te takes in the distance matrix and returns
  #                the standardized distance matrix. Defaults to the
  #                funcion that subtracks the minimum and divides by
  #                the range.
  #  caliper_type: Whether we want the caliper to be on DAPS or on the PS.
  #                caliper_type must either be 'DAPS', or 'PS'.
  #  coord_dist:   Set to true when we want to use a distance function that
  #                calculates the spherical distance of points instead of
  #                euclidean. Defaults to FALSE.
  #
  # Returns:
  #  List of next interval, matched dataset, standardized difference of
  #  the columns in cov.cols, indeces of matched treated and controls,
  #  and whether balance was achieved.
  
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




DAPSest <- function(dataset, out.col = NULL, trt.col = NULL, caliper = 0.1,
                    weight = 'optimal', coords.columns = NULL, ignore.cols = NULL,
                    pairsRet = FALSE, cov.cols = NULL, cutoff = 0.1, w_tol = 0.05,
                    distance = StandDist, caliper_type, quiet = FALSE,
                    coord_dist = FALSE, true_value = NULL) {
  
  # Function that fits DAPS by choosing the 'optimal' weight, or using
  # the weight specified and returns the causal effect estimate under a
  # specification of caliper. If the weight is set to 'optimal', the
  # optimal weight is chosen as the smallest weight that satisfies
  # standardized difference less than a cutoff. For computational reasons,
  # we do not fit DAPS on a lot of weights, but perform an iterative
  # algorithm to identify the optimal w.
  #
  # Args:
  #  dataset:       Data frame including treatment, outcome, coordinates,
  #                 and observed confounders.
  #  out.col:       If outcome column name is not 'Y', out.col should be
  #                 the index of the outcome column.
  #  trt.col:       If treatment is not named 'X', trt.col should be set
  #                 to the index of the treatment column.
  #  caliper:       A caliper for the DAPS Score difference of matched
  #                 pairs. Defaults to 0.1. Scalar.
  #  weight:        Scalar between 0 and 1 or should be set to 'optimal'.
  #                 Describes the percent of weight to be given on PS
  #                 difference. 1 - weight is given to standardized
  #                 distance. If set to 'optimal', the 'optimal' weight is
  #                 chosen. Defaults to 'optimal'.
  #  coords.columns:If the columns of coordinates are not named 'Longitude'
  #                 and 'Latitude', coords.columns are the column indeces
  #                 corresponding to longitude and latitude accordingly.
  #  ignore.cols:   All column indeces that should not be included in the
  #                 linear model. Often, this should be set to all columns
  #                 indeces corresponding to columns other than outcome,
  #                 treatment, and observed confounders.
  #  pairsRet:      Whether we want to return the information on the
  #                 matched pairs. Logical. Defaults to FALSE.
  #  cov.cols:      If the weight is set to 'optimal', standardized
  #                 difference of means will be calculated on the columns
  #                 whose indeces are in cov.cols.
  #  cutoff:        The cutoff of standardized difference of means under
  #                 which the covariates are considered balanced. Specify
  #                 when weight is set to 'optimal'. Defaults to 0.1.
  #  w_tol:         Tolerance on the choice of the optimal weight. Only
  #                 needed when weight is 'optimal'. Defaults to 0.05.
  #  distance:      Function te takes in the distance matrix and returns
  #                 the standardized distance matrix. Defaults to the
  #                 function that subtracks the minimum and divides by the
  #                 range.
  #  caliper_type:  Whether we want the caliper to be on DAPS or on the PS.
  #                 caliper_type must either be 'DAPS', or 'PS'.
  #  quiet:         Whether we want to print the choice of weight in DAPS
  #                 optimal. Defauls to TRUE.
  #  coord_dist:    Set to true when we want to use a distance function that
  #                 calculates the spherical distance of points instead of
  #                 euclidean. Defaults to FALSE.
  #  true_value:    Numeric. If provided, an indicator of whether the CI covers
  #                 the true value is returned.
  #
  # Returns:
  #  A list. The first element of the list is a vector of length 2. The
  #  first element of the vector is the causal effect estimated from a
  #  linear model on the matched pairs, adjusting for no observed
  #  confounding. The second element is the causal effect estimate from
  #  a linear model on the DAPS matched pairs that adjusts for all
  #  covariates that are not in ignore.cols. The second element of the
  #  list is a matrix of information on the matched pairs.
  
  # What we return.
  r <- NULL
  
  dataset <- as.data.frame(dataset)
  # Naming outcome and treatment as 'Y', 'X'.
  dataset <- FormDataset(dataset, ignore.cols = NULL,
                         out.col = out.col, trt.col = trt.col)
  # Fitting DAPS.
  
  if (is.numeric(weight)) {
    daps.out <- dist.ps(treated = dataset[dataset$X == 1, ],
                        control = dataset[dataset$X == 0, ],
                        caliper = caliper, weight = weight,
                        coords.columns = coords.columns,
                        distance = distance,
                        caliper_type = caliper_type,
                        coord_dist = coord_dist)
    # Getting the matched pairs.
    pairs.out        <- daps.out$match
    names(pairs.out) <- rownames(daps.out)
    pairs.out        <- na.omit(pairs.out)
    pairs.daps  <- dataset[c(as.numeric(names(pairs.out)),
                             as.numeric(pairs.out)), ]
    r$weight <- weight
    
  } else if (weight == 'optimal') {
    daps.opt <- DAPSopt(dataset, caliper = caliper,
                        coords.cols = coords.columns,
                        cov.cols = cov.cols, cutoff = cutoff,
                        w_tol = w_tol, distance = distance,
                        caliper_type = caliper_type, quiet = quiet,
                        coord_dist = coord_dist)
    pairs.daps <- daps.opt$pairs
    r$weight <- daps.opt$weight
    r$stand_diff <- daps.opt$stand_diff
    r$ind_trt <- daps.opt$ind_trt
    r$ind_cnt <- daps.opt$ind_cnt
  }
  
  # Dropping the columns that will not be included.
  # Outcome and treatment columns have already been renamed.
  pairs_small <- FormDataset(pairs.daps, ignore.cols = ignore.cols)
  est     <- numeric(2)
  lmod1 <- lm(Y ~ X, data = pairs_small)
  lmod2 <- lm(Y ~ X + ., data = pairs_small)
  est[1]  <- lmod1$coef[2]
  est[2]  <- lmod2$coef[2]
  
  r$est <- est
  
  if (!is.null(true_value)) {
    se_est <- c(summary(lmod1)$coef[2, 2], summary(lmod2)$coef[2, 2])
    r$cover <- (abs(true_value) - est < qnorm(0.975) * se_est)
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
