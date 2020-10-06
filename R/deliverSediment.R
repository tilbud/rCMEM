#' Calculate delivered sediment over one or more types of tides
#'
#' This function calculates the delivered sediment for a marsh elevation by summarising between one and three types of tidal cycles. 
#' @param z a numeric, elevation of the marsh
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
#' @param nFloods a numeric, the number of tidal flooding events per year
#' @param meanSeaLevel a numeric, Mean Sea Level
#' @param meanHighWater a numeric, Mean High Water level
#' @param meanHighHighWater a numeric, Mean Higher High Water level
#' @param meanHighHighWaterSpring a numeric, Mean Higher High Spring Tide Water level
#' @param meanLowWater a numeric, Mean Low Water level
#' @param meanLowLowWater a numeric, Mean Lower Low Water level
#' @param meanLowLowWaterSpring a numeric, Mean Low Lower Spring Tide Water level
#' @param captureRate a numeric, the number of times a water column will clear per tidal cycle
#' @param floodTime.fn a function, specify the method used to calculate flooding time per tidal cycle
#'
#' @return a numeric, the sediment delivered over the course of a year
#' @export
deliverSediment <- function(z, suspendedSediment, nFloods = 705.79,
                            meanSeaLevel, meanHighWater, meanHighHighWater=NA, meanHighHighWaterSpring=NA,
                            meanLowWater=meanSeaLevel-meanHighWater, 
                            meanLowLowWater=meanSeaLevel-meanHighHighWater, 
                            meanLowLowWaterSpring=meanSeaLevel-meanHighHighWaterSpring,
                            captureRate,
                            floodTime.fn = floodTimeLinear) {


  # If all three tidal datums are present
  if (all(! is.na(c(meanHighWater, meanHighHighWater, meanHighHighWaterSpring)))) {
    # Create a data frame operation so we can use tidy functions to speed this up
    tidalCycles <- data.frame(datumHigh = c(meanHighWater, meanHighHighWater, meanHighHighWaterSpring), 
                              datumLow = c(meanLowWater, meanLowLowWater, meanLowLowWaterSpring),
                              nTides = c(0.5, 0.46497542, 0.03502458)) 
  } else if (all(! is.na(c(meanHighWater, meanHighHighWater)))) {
    # If only MHW and MHHW are present
    tidalCycles <- data.frame(datumHigh = c(meanHighWater, meanHighHighWater), 
                              datumLow = c(meanLowWater, meanLowLowWater),
                              nTides = c(0.5, 0.5)) 
  } else {
    # If only MHW is present
    tidalCycles <- data.frame(datumHigh = c(meanHighWater), 
                              datumLow = c(meanLowWater),
                              nTides = c(1)) 
  }
  
  # Convert fraction of tides to counts of tides based on input
  tidalCycles$nTides <- tidalCycles$nTides * nFloods
  
  tidalCycles <- tidalCycles %>%
    # Set tidal propoerties to 0 if surface is above tidal range in each case
    dplyr::mutate(nTides = ifelse(z>datumHigh, 0, nTides), # number of tides
                  tidalHeight = ifelse(z>datumHigh, 0, (datumHigh-z)*0.5), # Tidal height relative to surface
                  floodTime = floodTime.fn(z=z, # Call flood time function.
                              datumHigh=datumHigh, 
                              datumLow=datumLow),
                  # Call available sediment function
                  fractionCaptured = ifelse(floodTime < 1/captureRate, # if the sediment column IS NOT able to clear
                                             # available suspendedSediment is total possible caputre 
                                             captureRate*floodTime, 
                                             # if the sediment column IS able to clear
                                             1),
                  availableSediment = suspendedSediment * nTides * tidalHeight, # Calculate available sediment as a cumulative block of water
                  deliveredSediment = availableSediment * fractionCaptured) # Calculated delivered sediment
  
  # Sum delivered sediment accross tidal cycles
  totalDeliveredSediment <- sum(tidalCycles$deliveredSediment) 
  
  return(totalDeliveredSediment)
}