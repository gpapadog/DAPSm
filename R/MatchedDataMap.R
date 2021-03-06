#' Plot matches.
#' 
#' Plots a map of the US of the matched pairs with color corresponding to
#' which unit is treated and which is the control and lines connecting
#' matched pairs.
#' 
#' @param x A data frame where the number of rows corresponds to the number
#' of matches. There are at least four columns in the data frame and correspond
#' to longitude and latitude information of the treated and the control unit.
#' Within the same row, we must have the coordinate information of treated and
#' control that are matched to each other.
#' @param trt_coords Indeces of the columns including the longitude and latitude
#' (with this order) of the treated unit.
#' @param con_coords Indeces of the columns including the longitude and latitude
#' (with this order) of the control unit.
#' @param plot_title Title of the plot.
#' @param point_data Whether we want to print points of the treated and control
#' observations. Treated units are green points and control units are red points.
#' 
#' @return Plot of matched pairs.
#' 
#' @export
#' @examples
#' data('toyData2')
#' toyData2$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#'                             data = toyData2)$fitted.values
#' daps1 <- DAPSest(toyData2, out.col = 2, trt.col = 1, caliper = 0.5,
#'                  weight = 0.2, coords.columns = c(4, 5),
#'                  pairsRet = TRUE, cov.cols = 6:9, cutoff = 0.1,
#'                  w_tol = 0.001, coord_dist = TRUE, caliper_type = 'DAPS',
#'                  matching_algorithm = 'greedy')
#' MatchedDataMap(x = daps1$pairs, trt_coords = c(3, 4), con_coords = c(7, 8))
#' # For a larger weight, less weight is given to distance, so matches will
#' # further apart.
#' daps2 <- DAPSest(toyData2, out.col = 2, trt.col = 1, caliper = 0.5,
#'                  weight = 0.8, coords.columns = c(4, 5),
#'                  pairsRet = TRUE, cov.cols = 6:9, cutoff = 0.1,
#'                  w_tol = 0.001, coord_dist = TRUE, caliper_type = 'DAPS')
#' MatchedDataMap(x = daps2$pairs, trt_coords = c(3, 4), con_coords = c(7, 8))
MatchedDataMap <- function(x, trt_coords, con_coords, plot.title = '',
                           point_data = TRUE) {
  
  x <- as.data.frame(x)
  
  point.data <- data.frame(lon = c(x[, trt_coords[1]], x[, con_coords[1]]),
                           lat = c(x[, trt_coords[2]], x[, con_coords[2]]))
  point.data$col <- 1
  point.data$col[1:(nrow(point.data) / 2)] <- 0
  point.data$col <- as.factor(point.data$col)
  
  line_data <- matrix(NA, nrow = 2 * nrow(x), ncol = 3)
  colnames(line_data) <- c('lon', 'lat', 'group')
  for (ii in 1:nrow(x)) {
    line_data[(2 * ii - 1) : (2 * ii), 1] <- as.numeric(x[ii, c(3, 7)])
    line_data[(2 * ii - 1) : (2 * ii), 2] <- as.numeric(x[ii, c(3, 7) + 1])
    line_data[(2 * ii - 1) : (2 * ii), 3] <- rep(ii, 2)
  }
  line_data <- as.data.frame(line_data)
  
  us.dat <- ggplot2::map_data("state")
  ct.dat <- ggplot2::map_data("county")
  
  g <- ggplot2::ggplot() +
    ggplot2::geom_polygon(ggplot2::aes(long, lat, group = group), color = 'grey55',
                          fill = 'grey85', data = us.dat) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.position = 'none',
                   panel.grid.minor = ggplot2::element_blank())
  if (point_data) {
    g <- g + ggplot2::geom_point(ggplot2::aes(x = lon, y = lat, color = col),
                                 data = point.data)
  }
  g <- g + ggplot2::scale_color_manual(values = c('darkgreen', 'red')) +
    ggplot2::geom_line(ggplot2::aes(x = lon, y = lat, group = group), data = line_data,
                       size = ggplot2::rel(0.4), col = 'black')
  if (plot.title != '') {
    g <- g + ggplot2::ggtitle(plot.title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = ggplot2::rel(1.6),
                                                        hjust = 0.5),
                     legend.key.size = unit(1, "cm"),
                     legend.title = ggplot2::element_text(size = ggplot2::rel(1)))
  }
  return(g)
}


