#' Calculate delivered sediment, three tidal cycle method
#'
#' This function calculates the delivered sediment for a marsh elevation over by summarising over three types of tidal cycles. 
#' @param z a numeric, elevation of the marsh
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
#' @param meanSeaLevel a numeric, Mean Sea Level
#' @param meanHighWater a numeric, Mean High Water level
#' @param meanHighHighWater a numeric, Mean Higher High Water level
#' @param meanHighHighWaterSpring a numeric, Mean Higher High Spring Tide Water level
#' @param meanLowWater a numeric, Mean Low Water level
#' @param meanLowLowWater a numeric, Mean Lower Low Water level
#' @param meanLowLowWaterSpring a numeric, Mean Low Lower Spring Tide Water level
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param capturedSediment a numeric, the fraction of available sediment captured by the marsh
#'
#' @return a numeric, the sediment delivered over the course of a year
#' @export
deliverSediment3TidalCycle <- function(z, suspendedSediment, meanSeaLevel, meanHighWater, meanHighHighWater, meanHighHighWaterSpring=NA,
                                         meanLowWater=meanSeaLevel-meanHighWater, meanLowLowWater=meanSeaLevel-meanHighHighWater, meanLowLowWaterSpring=meanSeaLevel-meanHighHighWaterSpring,
                                         settlingVelocity, capturedSediment=1) {
  
  if (! is.na(meanHighHighWaterSpring)) {
    # If MHHWS, MHHW, and MHW
    # Define Constants
    highTidesPerYear <- 352.657
    higherHighTidesPerYear <- 352.657 - 24.720
    springTidesPerYear <- 24.720
    hoursInTidalCycle <- 12.42
    
    # Create a data frame operation so we can use tidy functions to speed this up
    tidalCycles <- data.frame(datumHigh = c(meanHighWater, meanHighHighWater, meanHighHighWaterSpring), 
                              datumLow = c(meanLowWater, meanLowLowWater, meanLowLowWaterSpring),
                              nTides = c(highTidesPerYear, higherHighTidesPerYear, springTidesPerYear))
    
  } else {
    # If MHW and MHHW
    highTidesPerYear <- 352.657
    higherHighTidesPerYear <- 352.657
    hoursInTidalCycle <- 12.42
    
    # Create a data frame operation so we can use tidy functions to speed this up
    tidalCycles <- data.frame(datumHigh = c(meanHighWater, meanHighHighWater), 
                              datumLow = c(meanLowWater, meanLowLowWater),
                              nTides = c(highTidesPerYear, higherHighTidesPerYear))
    
  }
  
  tidalCycles <- tidalCycles %>%
    # Set tidal propoerties to 0 if surface is above tidal range in each case
    dplyr::mutate(nTides = ifelse(z>datumHigh, 0, nTides), # number of tides
                  tidalHeight = ifelse(z>datumHigh, 0, (datumHigh-z)*0.5), # Tidal height relative to surface
                  floodTime = ifelse(z>datumHigh, 0, floodTimeFromDatum(z=z, # Call flood time function.
                                                                        datumHigh=datumHigh, 
                                                                        datumLow=datumLow)), 
                  floodPct = ifelse(z>datumHigh, 0, floodTime / hoursInTidalCycle), # Convert hours to fraction
                  # Call available sediment function
                  availableSSC = availableSediment(floodPct=floodPct, suspendedSediment=suspendedSediment, 
                                                   settlingVelocity=settlingVelocity, capturedSediment=capturedSediment),
                  deliveredSSC = availableSSC * nTides * tidalHeight) # Calculated delivered SSC
  
  # Sum delivered SSC accross tidal cycles
  totalDeliveredSSC <- sum(tidalCycles$deliveredSSC) 
  
  return(totalDeliveredSSC)
}