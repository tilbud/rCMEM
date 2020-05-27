#' Calculate available sediment
#'
#' This function calculates the available sediment for a marsh elevation over a tidal cycle. 
#'
#' @param floodPct a numeric, fraction of time per tidal cycle the marsh is inundated
#' @param ssc a numeric, suspended sediment concentration of the water column
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param k numeric, the fraction of available sediment captured by the marsh
#'
#' @return numeric, the available sediment given a set of sediment and flooding conditions 
#' @export
availableSediment <- function(floodPct, ssc, settlingVelocity, k=1) {
  
  availableSSC <- ifelse(floodPct < 1/settlingVelocity, # if the sediment column IS NOT able to clear
                         # available ssc is total possible caputre 
                         ssc * floodPct/(1/settlingVelocity) * k, 
                         # if the sediment column IS able to clear
                         ssc * k)
  
  return(availableSSC)
}
