#' The central DAPS function
#'
#' Takes in a two data frames one with treatment and one with control units
#' including the variables: Longitude, Latitude and propensity scores (as
#' prop.scores) and returns a matrix of the matched pairs using DAPS, and
#' information on PS difference, matching difference, distance. Caliper can
#' can be specified on PS difference or DAPS.
#'
#' @param treated        A data frame include the treated units and the
#'                       variables: 'Longitude', 'Latitude' and propensity
#'                       scores (named 'prop.scores'). The rownames of treated
#'                       should be the unit ids.
#' @param control        Control units. Same variables as in treated.
#' @param caliper        A caliper of DAPS or PS difference for matching.
#' @param weight         Number between 0 and 1, percentage of matching weight
#'                       to be given on propensity score difference.
#' @param coords.columns If the columns of coordinates are not named
#'                      'Longitude', 'Latitude', coords.cols should be the
#'                       column indices corresponding to longitude and
#'                       latitude accordingly.
#' @param distance       Function that takes in the distance matrix and returns
#'                       the standardized distance matrix. Defaults to the
#'                       function that subtracks the minimum and divides by
#'                       the range.
#' @param caliper_type   Whether we want the caliper to be on DAPS or on the PS.
#'                       caliper_type must either be 'DAPS', or 'PS'.
#' @param coord_dist     Set to true when we want to use a distance function that
#'                       calculates the spherical distance of points instead of
#'                       euclidean. Defaults to FALSE.
#'
#' @return A dataframe, where each row corresponds to each treated unit, and
#' includes: the control unit to which it was matched, their propensity score
#' difference, their DAPS difference, their distance, their standardized distance.
#' NAs in the data frame correspond to units that were not matched.
#'
#' @examples
#' data('toyData')
#' toyData$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                            data = toyData)$fitted.values
#' daps <- dist.ps(treated = toyData[toyData$Z == 1, ],
#'                 control = toyData[toyData$Z == 0, ],
#'                 caliper_type = 'DAPS', caliper = 1,
#'                 coords.columns = c(4, 5))
#' head(daps)
dist.ps <- function(treated, control, caliper = 0.1, weight = 0.8,
                    coords.columns = NULL, distance = StandDist,
                    caliper_type = c('DAPS', 'PS'), coord_dist = FALSE) {
  
  caliper_type <- match.arg(caliper_type)

  if (!is.null(coords.columns)) {
    names(treated)[coords.columns] <- c('Longitude', 'Latitude')
    names(control)[coords.columns] <- c('Longitude', 'Latitude')
  }
  
  if (coord_dist) {  # Spherical distance on the globe.
    dist.mat <- fields::rdist.earth(cbind(treated$Longitude, treated$Latitude),
                                    cbind(control$Longitude, control$Latitude))
  } else {
    dist.mat <- fields::rdist(cbind(treated$Longitude, treated$Latitude),
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


