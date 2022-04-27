#' Calculate flood time using linear method
#'
#' This function calculates flood time in hours from an elevation, and a high and low tidal datum. 
#' @param z a numeric, elevation of the marsh
#' @param datumHigh a numeric, highest water level of the tidal period
#' @param datumLow a numeric, lowest water level of the tidal period
#' @param tidalCycleLength a numeric, the length of time of the tidal cycle. This is set to 1 to estimate flood time as a fraction of the tidal cycle.
#' 
#' @return a numeric, flood time per tidal period
#' @export
floodTimeLinear <- function(z, datumHigh, datumLow, tidalCycleLength = 1) {
  floodFract <- ifelse(z >= datumHigh, 0, 
                       ifelse(z <= datumLow, 1, 
                              (datumHigh-z)/(datumHigh-datumLow)))
  floodTime <- floodFract * tidalCycleLength
  return(floodTime)
}
