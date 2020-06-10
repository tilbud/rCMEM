#' Calculate delivered sediment, simple method
#'
#' This function calculates the delivered sediment for a marsh elevation over a year using the simplest possible approach. 
#' @param z a numeric, elevation of the marsh
#' @param floodPct a numeric, fraction of time per tidal cycle the marsh is inundated
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
#' @param meanHighWater a numeric, Mean High Water level
#' @param meanSeaLevel a numeric, Mean Sea Level
#' @param meanLowWater a numeric, Mean Low Water level
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param capturedSediment a numeric, the fraction of available sediment captured by the marsh
#' @param nTides an integer, the number of flooding events in a year
#'
#'
#' @return a numeric, the sediment delivered over the course of a year
#' @export 
deliveredSedimentSimple <- function(z, suspendedSediment, meanSeaLevel, meanHighWater, meanLowWater = meanSeaLevel-meanHighWater, settlingVelocity, 
                                    capturedSediment=1, nTides=704) {
  
  # Flood depth in meters is the same as water volume when in m and assuming a 1m2 area of interest
  floodDepth <- ifelse(z<=meanHighWater, (meanHighWater-z)*0.5, 0)
  
  # Mean flood percent is a line, relative position of elevation between meanHighWater and meanLowWater.
  floodPct <- ifelse(z >= meanHighWater, 0, ifelse(z <= meanLowWater, 1, (meanHighWater-z)/(meanHighWater-meanLowWater)))
  
  availableSSC <- availableSediment(suspendedSediment=suspendedSediment, floodPct=floodPct, settlingVelocity=settlingVelocity, capturedSediment=capturedSediment)
  
  # Mean flood depth is depth below meanHighWater. If above meanHighWater the marsh does not flood.
  # delivered sediment is the suspendedSediment available to the surface (g/m3) multiplied by 
  # the cumulative water volume (cube of water) and number of times the cube passes over the marsh. 
  deliveredSSC <- availableSSC * floodDepth * nTides
  
  return(deliveredSSC)
}

