#' Calculate flood time from datum
#'
#' This function calculates flood time in hours from an elevation, and a high and low tidal datum. 
#' @param z a numeric, elevation of the marsh
#' @param datumHigh a numeric, highest water level of the tidal period
#' @param datumLow a numeric, lowest water level of the tidal period
#' 
#' @references Hickey et al In Press. Tidal inundation modeling within GIS.
#' @references Department of the Navy. Hydrographic Branch. 1994. Australian National Tide Tables: Australia, Papua & New Guinea, Australian Hydrographic Publication 11.
#'  
#' @return a numeric, flood time in hours per tidal period
#' @export
floodTimeFromDatum <- function(z, datumHigh, datumLow) {
  
  # If elevation is above the tidal range indation time is 0
  datumHigh <- ifelse(z>=datumHigh, z, datumHigh)
  
  # If elevation is below inundation time is a full tidal cycle
  datumLow <- ifelse(z<=datumLow, z, datumLow)
  
  # Rising time over cell = 6.21 (A/pi - 1)	
  # where A = 2* pi - cos-1 [2 (height of cell – MLW) / (MHW – MLW) - 1] radians
  A1 <- 2 * pi - acos(2 * (z-datumLow) / (datumHigh-datumLow) - 1)
  risingTime <- 6.21 * (A1/pi - 1)
  
  # Falling time over cell = 6.21 (A/pi - 1) where
  # A = 2* - cos-1 [2 (height of cell – MHW) / (MLW – MHW) - 1] radians
  A2 <- 2 * pi - acos(2 * (z-datumHigh) / (datumLow-datumHigh) - 1)
  fallingTime <- 6.21 * (A2/pi - 1)
  
  # If between inundation time = abs (time rising - 6.21) + time falling
  inundationTime <- abs(risingTime - 6.21) + fallingTime
  
  return(inundationTime)
}
