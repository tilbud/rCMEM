#' Run Marsh Equilibrium Model with Cohorts Accross an Elevation Transect
#'
#' This function takes an elevation minimum, maximum, and intervals,  uses them to simulate an elevation transect, then batch runs runCohortMem over it.
#' @param startYear an integer, year in form YYYY, the start year of the scenario 
#' @param endYear an integer, year in form YYYY, the end year of the scenario  
#' @param relSeaLevelRiseInit a numeric, initial rate of relative sea-level rise
#' @param relSeaLevelRiseTotal a numeric, total relative sea-level rise over the course of the scanario
#' @param meanSeaLevel a numeric or a vector, Mean Sea Level at the start of the scenario, or a vector of Mean Sea Levels the same length as the number of years in a scenario
#' @param meanSeaLevelDatum a numeric, Mean Sea level over the last datum period
#' @param meanHighWaterDatum a numeric, Mean High Water level over the last datum period
#' @param meanHighHighWaterDatum a numeric (optional), Mean Higher High Water level over the last datum period
#' @param meanHighHighWaterSpringDatum a numeric (optional), Mean Higher High Spring Tide Water level over the last datum period
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
#' @param lunarNodalAmp a numeric, the amplitude of the 18-year lunar nodal cycle
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param bMax a numeric, or vector of numerics, maximum biomass
#' @param zVegMax a numeric, or vector of numerics, upper elevation of biomass limit
#' @param zVegMin a numeric, or vector of numerics, lower elevation of biomass limit
#' @param zVegPeak (optional) a numeric, or vector of numerics, elevation of peak biomass
#' @param plantElevationType character, either "orthometric" or "dimensionless", specifying elevation reference of the vegetation growing elevations
#' @param rootToShoot a numeric, or vector of numerics, root to shoot ratio
#' @param rootTurnover a numeric, or vector of numerics, belowground biomass annual turnover rate
#' @param abovegroundTurnover (optional) a numeric, or vector of numerics, aboveground biomass annual turnover rate
#' @param rootDepthMax a numeric, or vector of numerics, maximum (95\%) rooting depth
#' @param shape a character, "linear" or "exponential" describing the shape of the relationship between depth and root mass
#' @param omDecayRate a numeric, annual fractional mass lost
#' @param recalcitrantFrac a numeric, fraction of organic matter resistant to decay
#' @param omPackingDensity a numeric, the bulk density of pure organic matter
#' @param mineralPackingDensity a numeric, the bulk density of pure mineral matter
#' @param rootPackingDensity a numeric, the bulk density of pure root matter
#' @param speciesCode (optional) a character, or vector of characters, species names or codes associated with biological inputs
#' @param initialCohorts a data frame, (optional) custom set of mass cohorts that will override any decision making this function does
#' @param uplandCohorts a data frame,  (optional) custom set of mass cohorts to be used if intiial elevation is higher than both maximum tidal height and maximum wetland vegetation tolerance
#' @param supertidalCohorts a data frame, (optional) custom set of mass cohorts to be used if intiial elevation is higher than maximum tidal height, but not maximum wetland vegetation tolerance
#' @param supertidalSedimentInput a numeric, (optional) grams per cm^2 per year, an optional parameter which will define annual suspended sediment delivery to a sediment column that is is higher than maximum tidal height, but not maximum wetland vegetation tolerance
#' @param initElvMin a numeric, the lowest initial elevation over which to batch run the model
#' @param initElvMax a numeric, the highest initial elevation over which to batch run the model
#' @param elvIntervals a numeric, the elevation intervals of the transect
#' @param ...
#' 
#' @return a list of data frames, including the annualized summaries, and mapped cohorts tracked for every year, and every depth interval, of the simulation.
#' @export
runCohortMemTransect <- function(startYear, endYear=startYear+99, relSeaLevelRiseInit, relSeaLevelRiseTotal,
                                             meanSeaLevel, meanSeaLevelDatum=meanSeaLevel[1], 
                                             meanHighWaterDatum, meanHighHighWaterDatum=NA, meanHighHighWaterSpringDatum=NA, 
                                             suspendedSediment, lunarNodalAmp, settlingVelocity,
                                             bMax, zVegMin, zVegMax, zVegPeak, plantElevationType,
                                             rootToShoot, rootTurnover, rootDepthMax, shape="linear",
                                             omDecayRate, recalcitrantFrac,
                                             omPackingDensity=0.085, mineralPackingDensity=1.99,
                                             rootPackingDensity=omPackingDensity,
                                             abovegroundTurnover=NA,
                                             speciesCode=NA,
                                             initialCohorts=NA,
                                             uplandCohorts=NA,
                                             supertidalCohorts=NA,
                                             supertidalSedimentInput=NA,
                                             initElvMin=meanSeaLevel[1], 
                                             initElvMax=max(meanHighWaterDatum,meanHighHighWaterDatum, meanHighHighWaterSpringDatum, na.rm=T)+relSeaLevelRiseTotal, 
                                             elvIntervals=4, 
                                             ...) {
  
  elvTransect <- seq(initElvMin, initElvMax, elvIntervals)
  elvMins <- elvTransect - elvIntervals/2
  elvMaxs <- elvTransect + elvIntervals/2
  
  for (i in 1:length(elvTransect)) {
    memOutputs <- runCohortMem(initElv = elvTransect[i],
                                    startYear=startYear, 
                                    endYear=endYear, 
                                    relSeaLevelRiseInit=relSeaLevelRiseInit, 
                                    relSeaLevelRiseTotal=relSeaLevelRiseTotal,
                                    meanSeaLevel=meanSeaLevel, 
                                    meanSeaLevelDatum=meanSeaLevelDatum, 
                                    meanHighWaterDatum=meanHighWaterDatum, 
                                    meanHighHighWaterDatum=meanHighHighWaterDatum, 
                                    meanHighHighWaterSpringDatum=meanHighHighWaterSpringDatum, 
                                    suspendedSediment=suspendedSediment, 
                                    lunarNodalAmp=lunarNodalAmp,
                                    bMax=bMax, 
                                    zVegMin=zVegMin, 
                                    zVegMax=zVegMax, 
                                    zVegPeak=zVegPeak, 
                                    plantElevationType=plantElevationType,
                                    rootToShoot=rootToShoot, rootTurnover=rootTurnover, 
                                    rootDepthMax=rootDepthMax, shape=shape,
                                    omDecayRate=omDecayRate, recalcitrantFrac=recalcitrantFrac, 
                                    settlingVelocity=settlingVelocity,
                                    omPackingDensity=omPackingDensity, mineralPackingDensity=mineralPackingDensity,
                                    rootPackingDensity=rootPackingDensity,
                                    abovegroundTurnover=abovegroundTurnover,
                                    initialCohorts=initialCohorts,
                                    uplandCohorts=uplandCohorts,
                                    supertidalCohorts=supertidalCohorts,
                                    supertidalSedimentInput=supertidalSedimentInput)
    
    if (i == 1) {
      scenarioTransect <- memOutputs[[1]] %>% dplyr::mutate(elvMin = elvMins[i],
                                                     elvMax = elvMaxs[i],
                                                     initElv = elvTransect[i])
      cohortsTransect <-  memOutputs[[2]] %>% dplyr::mutate(elvMin = elvMins[i],
                                                     elvMax= elvMaxs[i],
                                                     initElv = elvTransect[i])
    } else {
      memOutputs[[1]] <- memOutputs[[1]] %>% dplyr::mutate(elvMin = elvMins[i],
                                 elvMax = elvMaxs[i],
                                 initElv = elvTransect[i])
      
      memOutputs[[2]] <- memOutputs[[2]] %>% dplyr::mutate(elvMin = elvMins[i],
                                 elvMax= elvMaxs[i],
                                 initElv = elvTransect[i])
      
      scenarioTransect <- scenarioTransect %>% 
        dplyr::bind_rows(memOutputs[[1]])
      
      cohortsTransect <- cohortsTransect %>% 
        dplyr::bind_rows(memOutputs[[2]])
    }
    
  }

  return(list(scenarioTransect=scenarioTransect,
              cohortsTransect=cohortsTransect))

}
