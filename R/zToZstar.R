#' Convert elevation to Z*
#'
#'  (z-MSL)/(MHW-MSL) 
#'
#' @param z a numeric elevation
#' @param MHW a numeric, mean high water tide
#' @param MSL a numeric, mean sea level
#'
#' @return elevation normalized to the tidal range
#' @export
#'
#' @examples
zToZstar <- function(z, MHW, MSL) { 
  
  (z-MSL)/(MHW-MSL) 
  
}