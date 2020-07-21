#' Convert elevation to Z*
#'
#'  (z-meanSeaLevel)/(meanHighWater-meanSeaLevel) 
#'
#' @param z a numeric elevation
#' @param meanHighWater a numeric, mean high water tide
#' @param meanSeaLevel a numeric, mean sea level
#'
#' @return elevation normalized to the tidal range
#' @export
#'
convertZToZstar <- function(z, meanHighWater, meanSeaLevel) { 
  
  (z-meanSeaLevel)/(meanHighWater-meanSeaLevel) 
  
}