#' Calculate flood time using trigonometric method
#'
#' This function calculates flood time in hours from an elevation, and a high and low tidal datum. 
#' @param z a numeric, elevation of the marsh
#' @param datumHigh a numeric, highest water level of the tidal period
#' @param datumLow a numeric, lowest water level of the tidal period
#' @param tidalCycleLength a numeric, the length of time of the tidal cycle. This is set to 1 so that it estimates fractional flood time
#' 
#' @references Hickey, R. (2019). Tidal inundation modeling within GIS. Journal of Coastal Conservation, 23(3), 599-606.
#' @references Department of the Navy. Hydrographic Branch. 1994. Australian National Tide Tables: Australia, Papua & New Guinea, Australian Hydrographic Publication 11.
#'  
#' @return a numeric, flood time per tidal period
#' @export
floodTimeTrig <- function(z, datumHigh, datumLow, tidalCycleLength = 1) {
  
  # If elevation is above the tidal range indation time is 0
  datumHigh <- ifelse(z>=datumHigh, z, datumHigh)
  
  # If elevation is below inundation time is a full tidal cycle
  datumLow <- ifelse(z<=datumLow, z, datumLow)
  
  # Rising time over cell = 6.21 (A/pi - 1)	
  # where A = 2* pi - cos-1 [2 (height of cell – meanLowWater) / (meanHighWater – meanLowWater) - 1] radians
  A1 <- 2 * pi - acos(2 * (z-datumLow) / (datumHigh-datumLow) - 1)
  risingTime <- tidalCycleLength/2 * (A1/pi - 1)
  
  # Falling time over cell = 6.21 (A/pi - 1) where
  # A = 2* - cos-1 [2 (height of cell – meanHighWater) / (meanLowWater – meanHighWater) - 1] radians
  A2 <- 2 * pi - acos(2 * (z-datumHigh) / (datumLow-datumHigh) - 1)
  fallingTime <- tidalCycleLength/2 * (A2/pi - 1)
  
  # If between inundation time = abs (time rising - 6.21) + time falling
  inundationTime <- abs(risingTime - tidalCycleLength/2) + fallingTime
  
  return(inundationTime)
}
