#' Annual mineral inputs into marsh surface
#' 
#' This function cacluates the annaual sediment input into the marsh surface via the following formula:
#' Number of tides per year, times the unit area of interest, times the mean tidal volume above the marsh surface in cm$^3$, times the suspended sediment concentration in mg per liter, times 1e-6 conversion factor to convert from liter to cm$^3, resulting in the mineral mass added per year. The number of times (nTidesPerYear), unit area length (soilLength), and unit area width (soilWidth) are included in the consts list (see specific names).
#'
#' @param sse a numeric representing suspended sediment concentration in (mg per liter)
#' @param meanTidalVolume a numeric representing the the mean tide volume above the marsh surface (cm^33)
#' @param parms placeholder, currently NULL
#' @param consts a list of constants including nTidesPerYear, soilLength, and soilWidth
#'
#' @return a numeric that is the mass of mineral added in one year to the top of the marsh
#' @export
#' 
#' @example 
#' sedimentInputs(ssc=20, meanTidalHeight=10, consts=list(nTidesPerYear = 704, soilLength = 1, soilWidth = 1))
#'
sedimentInputs <- function(ssc, # Suspended Sediment Concentration, mg per liter
                           meanTidalHeight, #mean tide height above marsh
                           parms=NULL, 
                           consts){ 
  ##Check constants
  if(!all(c('nTidesPerYear', 'soilLength', 'soilWidth') %in% names(consts))){
    stop('Can not find expected constant names')
  }
  
  if(meanTidalHeight < 0){
    return(0) #only add sediment to the marsh if the mean high water is above the marsh elevation
  }
  
  # convert mg/l to grams/cm^3
  ssc_gPerCm3 <- ssc * 1e-6 
  
  # Cumulative water volume
  cumAnnWaterVol <- consts$nTidesPerYear * meanTidalHeight * consts$soilLength * consts$soilWidth
  
  annSediment <- ssc_gPerCm3 * cumAnnWaterVol #g-sediment per year
  
  return(annSediment)
}
