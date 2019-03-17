#' Annual mineral inputs into marsh
#'
#' @param sse a numeric representing suspended sediment concentration, generally mg per liter
#' @param meanTidalVolume a numeric representing the the mean tide volume above the marsh (cm3)
#' @param parms placeholder, currently NULL
#' @param consts a list of constants including nTidesPerYear
#'
#' @return a numeric that is the mass of mineral added in one year to the top of the marsh
#' @export
#'
#' @examples
sedimentInputs <- function(ssc, # Suspended Sediment Concentration, mg per liter
                           meanTidalHeight, #mean tide height above marsh
                           parms=NULL, consts){ 
  ##Check constants
  if(!all(c('nTidesPerYear', 'soilLength', 'soilWidth') %in% names(consts))){
    stop('Can not find expected constant names')
  }
  
  if(meanTidalHeight < 0){
    return(0) #only add sediment to the marsh if the mean high water is above the marsh elevation
  }
  
  ssc_gPerCm3 <- ssc * 1e-6 # convert mg/l to grams/cm^3
  
   # Cumulative water volume
  cumAnnWaterVol <- consts$nTidesPerYear * meanTidalHeight * consts$soilLength * consts$soilWidth
  
  annSediment <- ssc_gPerCm3 * cumAnnWaterVol #g-sediment per year
  
  return(annSediment)
}
