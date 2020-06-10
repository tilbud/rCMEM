#' Calculate available sediment
#'
#' This function calculates the available sediment for a marsh elevation over a tidal cycle. 
#'
#' @param floodPct a numeric, fraction of time per tidal cycle the marsh is inundated
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param capturedSediment numeric, the fraction of available sediment captured by the marsh
#'
#' @return numeric, the available sediment given a set of sediment and flooding conditions 
#' @export
availableSediment <- function(floodPct, suspendedSediment, settlingVelocity, capturedSediment=1) {
  
  availableSSC <- ifelse(floodPct < 1/settlingVelocity, # if the sediment column IS NOT able to clear
                         # available suspendedSediment is total possible caputre 
                         suspendedSediment * floodPct/(1/settlingVelocity) * capturedSediment, 
                         # if the sediment column IS able to clear
                         suspendedSediment * capturedSediment)
  
  return(availableSSC)
}
