#' Run Marsh Equilibrium Model with Cohorts
#'
#' This function takes an initial elevation, average annual suspended sediment concentrtion and tidal properties, plant traits, and decay rates, then builds a scenario, and runs the marsh equilibrium model over that scenario tracking individual sediment mass cohorts including minderal, live root, slow decay organic mater pool and a fast organic matter pool.
#' @param startYear an integer, year in form YYYY, the start year of the scenario 
#' @param endYear an integer, year in form YYYY, the end year of the scenario  
#' @param relSeaLevelRiseInit a numeric, initial rate of relative sea-level rise
#' @param relSeaLevelRiseTotal a numeric, total relative sea-level rise over the course of the scanario
#' @param initElv a numeric, the initial elevation of the marsh at the start of the scenario
#' @param meanSeaLevel a numeric or a vector, Mean Sea Level at the start of the scenario, or a vector of Mean Sea Levels the same length as the number of years in a scenario
#' @param meanSeaLevelDatum a numeric, Mean Sea level over the last datum period
#' @param meanHighWaterDatum a numeric, Mean High Water level over the last datum period
#' @param meanHighHighWaterDatum a numeric (optional), Mean Higher High Water level over the last datum period
#' @param meanHighHighWaterSpringDatum a numeric (optional), Mean Higher High Spring Tide Water level over the last datum period
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
#' @param lunarNodalAmp a numeric, the amplitude of the 18-year lunar nodal cycle
#' @param lunarNodalPhase a numeric, in decimal years (YYYY) the start year of the sine wave representing the lunar nodal cycle 
#' @param captureRate a numeric, the number of times a water column will clear per tidal cycle
#' @param nFloods a numeric, the number of tidal flooding events per year
#' @param floodTime.fn a function, specify the method used to calculate flooding time per tidal cycle
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
#' @param abovegroundTurnover (optional) a numeric or vector of numerics, aboveground biomass annual turnover rate
#' @param speciesCode (optional) a character, or vector of characters, species names or codes associated with biological inputs
#' @param initialCohorts a data frame, (optional) custom set of mass cohorts that will override any decision making this function does
#' @param uplandCohorts a data frame,  (optional) custom set of mass cohorts to be used if intiial elevation is higher than both maximum tidal height and maximum wetland vegetation tolerance
#' @param superTidalCohorts a data frame, (optional) custom set of mass cohorts to be used if intiial elevation is higher than maximum tidal height, but not maximum wetland vegetation tolerance
#' @param supertidalSedimentInput, a numeric, (optional) grams per cm^2 per year, an optional parameter which will define annual suspended sediment delivery to a sediment column that is is higher than maximum tidal height, but not maximum wetland vegetation tolerance
#' @param ...
#' 
#' @return a list of data frames, including the annualized summaries, and mapped cohorts tracked for every year of the simulation.
#' @export
runCohortMem <- function(startYear, endYear=startYear+99, relSeaLevelRiseInit, relSeaLevelRiseTotal, initElv,
                         meanSeaLevel, meanSeaLevelDatum=meanSeaLevel[1], meanHighWaterDatum, meanHighHighWaterDatum=NA, meanHighHighWaterSpringDatum=NA, 
                         suspendedSediment, lunarNodalAmp, lunarNodalPhase=2011.181,
                         nFloods = 705.79, floodTime.fn = floodTimeLinear,
                         bMax, zVegMin, zVegMax, zVegPeak, plantElevationType,
                         rootToShoot, rootTurnover, abovegroundTurnover=NA, speciesCode=NA, rootDepthMax, shape="linear",
                         omDecayRate, recalcitrantFrac, captureRate,
                         omPackingDensity=0.085, mineralPackingDensity=1.99,
                         rootPackingDensity=omPackingDensity,
                         initialCohorts=NA,
                         uplandCohorts=NA,
                         supertidalCohorts=NA,
                         supertidalSedimentInput=NA,
                         ...) {
  
  # Make sure tidyverse is there
  ##TODO move this to the package description
  require(tidyverse, quietly = TRUE)
  
  # Build scenario curve
  scenario <- buildScenarioCurve(startYear=startYear, endYear=endYear, 
                                 meanSeaLevel=meanSeaLevel, 
                                 relSeaLevelRiseInit=relSeaLevelRiseInit, 
                                 relSeaLevelRiseTotal=relSeaLevelRiseTotal, 
                                 suspendedSediment=suspendedSediment)
  
  # add high tides
  scenario <- buildHighTideScenario(scenario, 
                                    meanSeaLevelDatum=meanSeaLevelDatum, 
                                    meanHighWaterDatum=meanHighWaterDatum, 
                                    meanHighHighWaterDatum=meanHighHighWaterDatum, 
                                    meanHighHighWaterSpringDatum=meanHighHighWaterSpringDatum, 
                                    lunarNodalAmp = lunarNodalAmp,
                                    lunarNodalPhase=lunarNodalPhase)
  
  # Add blank colums for attributes we will add later
  scenario$surfaceElevation <- as.numeric(rep(NA, nrow(scenario)))
  scenario$speciesCode <- as.character(rep(NA, nrow(scenario)))
  scenario$rootToShoot <- as.numeric(rep(NA, nrow(scenario)))
  scenario$rootTurnover <- as.numeric(rep(NA, nrow(scenario)))
  scenario$abovegroundTurnover <- as.numeric(rep(NA, nrow(scenario)))
  scenario$rootDepthMax <- as.numeric(rep(NA, nrow(scenario)))
  scenario$aboveground_biomass <- as.numeric(rep(NA, nrow(scenario)))
  scenario$belowground_biomass <- as.numeric(rep(NA, nrow(scenario)))
  scenario$mineral <- as.numeric(rep(NA, nrow(scenario)))

  # Set initial conditions
  initialConditions <- determineInitialCohorts(initElv=initElv,
                                               meanSeaLevel=scenario$meanSeaLevel[1], 
                                               meanHighWater=scenario$meanHighWater[1], 
                                               meanHighHighWater=scenario$meanHighHighWater[1], 
                                               meanHighHighWaterSpring=scenario$meanHighHighWaterSpring[1],
                                               suspendedSediment=scenario$suspendedSediment[1],
                                               nFloods = nFloods, floodTime.fn = floodTime.fn,
                                               bMax=bMax, zVegMin=zVegMin, zVegMax=zVegMax, zVegPeak=zVegPeak, 
                                               plantElevationType=plantElevationType,
                                               rootToShoot=rootToShoot, rootTurnover=rootTurnover, 
                                               rootDepthMax=rootDepthMax, shape=shape,
                                               abovegroundTurnover=abovegroundTurnover,
                                               omDecayRate=omDecayRate, 
                                               recalcitrantFrac=recalcitrantFrac, 
                                               captureRate=captureRate,
                                               omPackingDensity=omPackingDensity, 
                                               mineralPackingDensity=mineralPackingDensity,
                                               rootPackingDensity=omPackingDensity,
                                               speciesCode=speciesCode,
                                               initialCohorts=initialCohorts,
                                               uplandCohorts=uplandCohorts,
                                               supertidalCohorts=supertidalCohorts,
                                               supertidalSedimentInput=supertidalSedimentInput)
  cohorts <- initialConditions[[1]]

  # Add initial conditions to annual time step tracker
  scenario$surfaceElevation[1] <- initElv
  scenario[1,names(scenario) %in% names(initialConditions[[2]])] <- initialConditions[[2]]
  scenario$mineral[1] <- initialConditions[[3]]
  
  cohorts$year <- rep(scenario$year[1], nrow(cohorts))
  
  # Preallocate memory for cohort tracking
  nInitialCohorts <- nrow(cohorts)
  nScenarioYears <- nrow(scenario)
  initCohortRows <- nInitialCohorts * nScenarioYears
  newCohortRows <- sum(1:nScenarioYears)
  totalRows <- initCohortRows+newCohortRows
  
  trackCohorts <- data.frame(age=as.numeric(rep(NA, totalRows)),
                             fast_OM=as.numeric(rep(NA, totalRows)),
                             slow_OM=as.numeric(rep(NA, totalRows)),
                             respired_OM=as.numeric(rep(NA, totalRows)),
                             mineral=as.numeric(rep(NA, totalRows)), 
                             root_mass=as.numeric(rep(NA, totalRows)),
                             layer_top=as.numeric(rep(NA, totalRows)),
                             layer_bottom=as.numeric(rep(NA, totalRows)),
                             cumCohortVol=as.numeric(rep(NA, totalRows)),
                             year=as.integer(rep(NA, totalRows)))
  
  # add initial set of cohorts
  trackCohorts[1:nInitialCohorts,] <- cohorts
  
  # create variables to keep track of cohorts added to the full cohort tracker
  cohortsNewRowMin <- nInitialCohorts + 1
  
  # Calculate the unmoving bottom of the profile as a consistent reference point
  profileBottomElv <- initElv - max(cohorts$layer_bottom)
  
  # Convert real growing elevations to dimensionless growing elevations
  if (! plantElevationType %in% c("dimensionless", "zStar", "Z*", "zstar")) {
    zStarVegMin <- convertZToZstar(zVegMin, meanHighWaterDatum, meanSeaLevelDatum)
    zStarVegMax <- convertZToZstar(zVegMax, meanHighWaterDatum, meanSeaLevelDatum)
    zStarVegPeak <- convertZToZstar(zVegPeak, meanHighWaterDatum, meanSeaLevelDatum)
  } else {
    zStarVegMin <- zVegMin
    zStarVegMax <- zVegMax
    zStarVegPeak <- zVegPeak
  }
  
  # Fourth, add one cohort for each year in the scenario
  # Iterate through scenario table
  for (i in 2:nrow(scenario)) {
    
    # Calculate surface elevation relative to datum
    surfaceElvZStar <- convertZToZstar(z=scenario$surfaceElevation[i-1], 
                                       meanHighWater=scenario$meanHighWater[i], 
                                       meanSeaLevel=scenario$meanSeaLevel[i])
    
    # Calculate dynamic above ground biomass
    bio_table <- runMultiSpeciesBiomass(z=surfaceElvZStar, bMax=bMax, zVegMax=zStarVegMax, 
                                        zVegMin=zStarVegMin, zVegPeak=zStarVegPeak,
                                        rootToShoot=rootToShoot,
                                        rootTurnover=rootTurnover, 
                                        abovegroundTurnover=abovegroundTurnover, 
                                        rootDepthMax=rootDepthMax, 
                                        speciesCode=speciesCode)    
    
    # Calculate Mineral pool
    dynamicMineralPool <- deliverSediment(z=scenario$surfaceElevation[i-1], 
                                          suspendedSediment=scenario$suspendedSediment[i], 
                                          meanSeaLevel=scenario$meanSeaLevel[i], 
                                          meanHighWater=scenario$meanHighWater[i], 
                                          meanHighHighWater = scenario$meanHighHighWater[i], 
                                          meanHighHighWaterSpring = scenario$meanHighHighWaterSpring[i], 
                                          captureRate=captureRate,
                                          nFloods=nFloods,
                                          floodTime.fn=floodTime.fn)
    
   
    # Add a the new inorganic sediment cohort,
    # add live roots, and age the organic matter
    cohorts <- addCohort(cohorts, totalRootMassPerArea=bio_table$belowground_biomass[1], rootDepthMax=bio_table$rootDepthMax[1], 
                         rootTurnover = bio_table$rootTurnover[1], omDecayRate = list(fast=omDecayRate, slow=0),
                         rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                         packing=list(organic=omPackingDensity, mineral=mineralPackingDensity), 
                         rootDensity=rootPackingDensity, shape=shape,
                         mineralInput = dynamicMineralPool)
    
    # Add the calendar year to the cohort table so that we can track each cohort over each year
    cohorts$year <- rep(scenario$year[i], nrow(cohorts))
    
    # add new set of cohorts, to the table of all the cohorts we are tracking over each year
    cohortsNewRowMax <- cohortsNewRowMin + nrow(cohorts) - 1
    trackCohorts[cohortsNewRowMin:cohortsNewRowMax,] <- cohorts
    
    # Keep track of the number of rows in the cohorts table
    cohortsNewRowMin <- cohortsNewRowMin + nrow(cohorts) + 1
    
    # Add annual variables to annual time step summary table
    scenario$surfaceElevation[i] <- profileBottomElv + max(cohorts$layer_bottom, na.rm=T)
    scenario[i,names(scenario) %in% names(bio_table)] <- bio_table
    scenario$mineral[i] <- dynamicMineralPool
    
  }
  
  # Remove NA values from cohorts
  trackCohorts <- trimCohorts(trackCohorts)
  
  # Return annual time steps and full set of cohorts
  outputsList <- list(annualTimeSteps = scenario, cohorts = trackCohorts)
  
  return(outputsList)
}