#' Calculate delivered sediment, simple method
#'
#' This function calculates the delivered sediment for a marsh elevation over a year using the simplest possible approach. 
#' @param z a numeric, elevation of the marsh
#' @param floodPct a numeric, fraction of time per tidal cycle the marsh is inundated
#' @param ssc a numeric, suspended sediment concentration of the water column
#' @param MHW a numeric, Mean High Water level
#' @param MSL a numeric, Mean Sea Level
#' @param MLW a numeric, Mean Low Water level
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param k a numeric, the fraction of available sediment captured by the marsh
#' @param nTides an integer, the number of flooding events in a year
#'
#'
#' @return a numeric, the sediment delivered over the course of a year
#' @export 
deliveredSedimentSimple <- function(z, ssc, MSL, MHW, MLW = MSL-MHW, settlingVelocity, 
                                    k=1, nTides=704) {
  
  # Flood depth in meters is the same as water volume when in m and assuming a 1m2 area of interest
  floodDepth <- ifelse(z<=MHW, (MHW-z)*0.5, 0)
  
  # Mean flood percent is a line, relative position of elevation between MHW and MLW.
  floodPct <- ifelse(z >= MHW, 0, ifelse(z <= MLW, 1, (MHW-z)/(MHW-MLW)))
  
  availableSSC <- availableSediment(ssc=ssc, floodPct=floodPct, settlingVelocity=settlingVelocity, k=k)
  
  # Mean flood depth is depth below MHW. If above MHW the marsh does not flood.
  # delivered sediment is the ssc available to the surface (g/m3) multiplied by 
  # the cumulative water volume (cube of water) and number of times the cube passes over the marsh. 
  deliveredSSC <- availableSSC * floodDepth * nTides
  
  return(deliveredSSC)
}

