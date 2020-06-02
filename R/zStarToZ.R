#' Convert from Z* to elevation
#'
#' @param zStar a numeric describing the normalized Z* elevation
#' @param meanHighWater a numeric, mean high water tide
#' @param meanSeaLevel a numeric, mean sea level
#'
#' @return a numerical describing the marsh elevation
#' @export
#'
zStarToZ <- function(zStar, meanHighWater, meanSeaLevel) { 
  (zStar * ((meanHighWater-meanSeaLevel))) + meanSeaLevel 
  }