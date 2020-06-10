#' Run Marsh Equilibrium Model with Cohorts Accross an Elevation Transect
#'
#' This function takes an initial elevation, average annual suspended sediment concentrtion and tidal properties, plant traits, and decay rates, then builds a scenario, and runs the marsh equilibrium model over that scenario tracking individual sediment mass cohorts including minderal, live root, slow decay organic mater pool and a fast organic matter pool.
#' @param initElvMin
#' @param initElvMax
#' @param elvIntervals
#' @param startYear an integer, year in form YYYY, the start year of the scenario 
#' @param endYear an integer, year in form YYYY, the end year of the scenario  
#' @param rslrT1 a numeric, initial rate of relative sea-level rise
#' @param rslrTotal a numeric, total relative sea-level rise over the course of the scanario
#' @param MSL a numeric or a vector, Mean Sea Level at the start of the scenario, or a vector of Mean Sea Levels the same length as the number of years in a scenario
#' @param MSL0 a numeric, Mean Sea level over the last datum period
#' @param MHW a numeric, Mean High Water level over the last datum period
#' @param MHHW a numeric (optional), Mean Higher High Water level over the last datum period
#' @param MHHWS a numeric (optional), Mean Higher High Spring Tide Water level over the last datum period
#' @param ssc a numeric, suspended sediment concentration of the water column
#' @param lunarNodalAmp a numeric, the amplitude of the 18-year lunar nodal cycle
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param bMax a numeric, maximum biomass
#' @param zVegMax a numeric, upper elevation of biomass limit
#' @param zVegMin a numeric, lower elevation of biomass limit
#' @param zVegPeak a numeric, (optional) elevation of peak biomass
#' @param plantElevationType character, either "orthometric" or "dimensionless", specifying elevation reference of the vegetation growing elevations
#' @param rootToShoot a numeric, root to shoot ratio
#' @param rootTurnover a numeric, below ground biomass annual turnover rate
#' @param rootDepthMax a numeric, maximum (95\%) rooting depth
#' @param shape a character, "linear" or "exponential" describing the shape of the relationship between depth and root mass
#' @param omDecayRate a numeric, annual fractional mass lost
#' @param recalcitrantFrac a numeric, fraction of organic matter resistant to decay
#' @param omPackingDensity a numeric, the bulk density of pure organic matter
#' @param mineralPackingDensity a numeric, the bulk density of pure mineral matter
#' @param rootPackingDensity a numeric, the bulk density of pure root matter
#' @param initialCohorts a data frame, (optional) custom set of mass cohorts that will override any decision making this function does
#' @param uplandCohorts a data frame,  (optional) custom set of mass cohorts to be used if intiial elevation is higher than both maximum tidal height and maximum wetland vegetation tolerance
#' @param superTidalCohorts a data frame, (optional) custom set of mass cohorts to be used if intiial elevation is higher than maximum tidal height, but not maximum wetland vegetation tolerance
#' @param supertidalSedimentInput, a numeric, (optional) grams per cm^2 per year, an optional parameter which will define annual suspended sediment delivery to a sediment column that is is higher than maximum tidal height, but not maximum wetland vegetation tolerance
#' @param ...
#' 
#' @return a list of data frames, including the annualized summaries, and mapped cohorts tracked for every year, and every depth interval, of the simulation.
#' @export
runMemWithCohortsAccrossTransect <- function(initElvMin=MSL[1], initElvMax=max(MHW,MHHW, MHHWS, na.rm=T)+rslrTotal, elvIntervals=4, 
                                             startYear, endYear=startYear+99, rslrT1, rslrTotal,
                                             MSL, MSL0=MSL[1], MHW, MHHW=NA, MHHWS=NA, ssc, lunarNodalAmp,
                                             bMax, zVegMin, zVegMax, zVegPeak, plantElevationType,
                                             rootToShoot, rootTurnover, rootDepthMax, shape="linear",
                                             omDecayRate, recalcitrantFrac, settlingVelocity,
                                             omPackingDensity=0.085, mineralPackingDensity=1.99,
                                             rootPackingDensity=omPackingDensity,
                                             initialCohorts=NA,
                                             uplandCohorts=NA,
                                             supertidalCohorts=NA,
                                             supertidalSedimentInput=NA,
                                             ...) {
  
  elvTransect <- seq(initElvMin, initElvMax, elvIntervals)
  elvMins <- elvTransect - elvIntervals/2
  elvMaxs <- elvTransect + elvIntervals/2
  
  for (i in 1:length(elvTransect)) {
    memOutputs <- runMemWithCohorts(initElv = elvTransect[i],
                      startYear=startYear, endYear=endYear, rslrT1=rslrT1, rslrTotal=rslrTotal,
                      MSL=MSL, MSL0=MSL0, MHW=MHW, MHHW=MHHW, MHHWS=MHHWS, ssc=ssc, lunarNodalAmp=lunarNodalAmp,
                      bMax=bMax, zVegMin=zVegMin, zVegMax=zVegMax, zVegPeak=zVegPeak, plantElevationType=plantElevationType,
                      rootToShoot=rootToShoot, rootTurnover=rootTurnover, rootDepthMax=rootDepthMax, shape=shape,
                      omDecayRate=omDecayRate, recalcitrantFrac=recalcitrantFrac, settlingVelocity=settlingVelocity,
                      omPackingDensity=omPackingDensity, mineralPackingDensity=mineralPackingDensity,
                      rootPackingDensity=rootPackingDensity,
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
  
  scenarioTransectGraph <- scenarioTransect %>% 
    filter(surfaceElevation != initElv) %>% 
    mutate(above_below_msl = ifelse(surfaceElevation >= MSL, "above", "below"))
  
  ggplot(data=scenarioTransectGraph, aes(x=years, y=surfaceElevation, color=above_below_msl)) +
    geom_line(aes(group=as.character(initElv))) +
    geom_point()
  
  return(list(scenarioTransect=scenarioTransect,
              cohortsTransect=cohortsTransect))

}
