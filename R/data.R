#' Toy data set.
#' 
#' A data set containing simulated data with a binary treatment, continuous
#' outcome, continuous covariates, and coordinates of the points.
#' 
#'
#' @format A data frame with 800 rows and 9 columns.
#' \describe{
#'   \item{Z}{A binary treatment. Values 0, 1.}
#'   \item{Y}{The continuous outcome.}
#'   \item{U}{Spatial covariate with Matern correlation function.}
#'   \item{long}{Longitude of the observation.}
#'   \item{lat}{Latitude of the observation.}
#'   \item{X1}{Continuous covariate simulated independently as N(0, 1).}
#'   \item{X2}{Continuous covariate simulated independently as N(0, 1).}
#'   \item{X3}{Continuous covariate simulated independently as N(0, 1).}
#'   \item{X4}{Continuous covariate simulated independently as N(0, 1).}
#' }
"toyData"


#' Toy data set.
#' 
#' A data set containing simulated data with a binary treatment, continuous
#' outcome, continuous covariates, and coordinates of the points. The locations
#' of the observations are power plant facilities.
#'
#' @format A data frame with 200 rows and 9 columns.
#' \describe{
#'   \item{Z}{A binary treatment. Values 0, 1.}
#'   \item{Y}{The continuous outcome.}
#'   \item{U}{Spatial covariate with Matern correlation function.}
#'   \item{long}{Longitude of the observation.}
#'   \item{lat}{Latitude of the observation.}
#'   \item{X1}{Continuous covariate simulated independently as N(0, 1).}
#'   \item{X2}{Continuous covariate simulated independently as N(0, 1).}
#'   \item{X3}{Continuous covariate simulated independently as N(0, 1).}
#'   \item{X4}{Continuous covariate simulated independently as N(0, 1).}
#' }
"toyData2"