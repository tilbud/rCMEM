#' Convert from Z* to elevation
#'
#' @param zStar a numeric describing the normalized Z* elevation
#' @param MHW a numeric, mean high water tide
#' @param MSL a numeric, mean sea level
#'
#' @return a numerical describing the marsh elevation
#' @export
#'
#' @examples
zStarToZ <- function(zStar, MHW, MSL) { 
  (zStar * ((MHW-MSL))) + MSL 
  }