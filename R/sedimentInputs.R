#' Annual mineral inputs into marsh surface
#' 
#' This function cacluates the annaual sediment input into the marsh surface via the following formula:
#' Number of tides per year, times the unit area of interest, times the mean tidal volume above the marsh surface in cm$^3$, times the suspended sediment concentration in mg per liter, times 1e-6 conversion factor to convert from liter to cm$^3, resulting in the mineral mass added per year. The number of times (nTidesPerYear), unit area length (soilLength), and unit area width (soilWidth) are included in the consts list (see specific names).
#'
#' @param meanTidalHeight annual mean tide above marsh elevation
#' @param nTidesPerYear number of tides per year
#' @param soilLength unit length of interest
#' @param soilWidth unit width of interest
#' @param ... 
#' @param ssc a numeric representing suspended sediment concentration in (mg per liter)
#'
#' @return a numeric that is the mass of mineral added in one year to the top of the marsh
#' @export
#' 
#'
sedimentInputs <- function(ssc, # Suspended Sediment Concentration, mg per liter
                           meanTidalHeight, #mean tide height above marsh
                           nTidesPerYear = 704,
                           soilLength=1, soilWidth=1, ...){ 
  if(meanTidalHeight < 0){
    return(0) #only add sediment to the marsh if the mean high water is above the marsh elevation
  }

  annSediment <- (ssc * 1e-6)  * # convert mg/l to grams/cm^3
    nTidesPerYear * (meanTidalHeight * 0.5 * soilLength * soilWidth) # Cumulative water volume
  
  if (annSediment < 0) {
    return(0)
  } else {
    return(annSediment)
  }
}
