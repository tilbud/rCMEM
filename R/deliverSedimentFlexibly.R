#' Deliver sediment flexibly
#'
#' Function decides which sediment module to use based on available inputs
#' @param z a numeric, elevation of the marsh
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
#' @param meanSeaLevel a numeric, Mean Sea Level
#' @param meanHighWater a numeric, Mean High Water level
#' @param meanHighHighWater a numeric, Mean Higher High Water level
#' @param meanHighHighWaterSpring a numeric, Mean Higher High Spring Tide Water level
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param ...
#' 
#' @return a numeric, sediment delivered in a year of flooding
#' @export
deliverSedimentFlexibly <- function(z, suspendedSediment, meanSeaLevel, meanHighWater, meanHighHighWater=NA, meanHighHighWaterSpring=NA, settlingVelocity, ...) {
  # If the user does not include detailed tidal datums
  deliveredSSC <- ifelse(all(is.na(c(meanHighHighWater, meanHighHighWaterSpring))),
                         # run the simple SSC module
                         deliverSedimentSimple(z=z, 
                                                 suspendedSediment=suspendedSediment, 
                                                 meanSeaLevel=meanSeaLevel, meanHighWater=meanHighWater, 
                                                 settlingVelocity=settlingVelocity),
                         # If they do include them, run the 3 tide stage module
                         deliverSediment3TidalCycle(z=z, 
                                                      suspendedSediment=suspendedSediment, 
                                                      meanSeaLevel=meanSeaLevel, 
                                                      meanHighWater=meanHighWater, 
                                                      settlingVelocity=settlingVelocity,
                                                      meanHighHighWater=meanHighHighWater, 
                                                      meanHighHighWaterSpring=meanHighHighWaterSpring))
  
  return(deliveredSSC)
  
}