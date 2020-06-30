#' Run Marsh Equilibrium Model with Cohorts
#'
#' This function takes an initial elevation, average annual suspended sediment concentrtion and tidal properties, plant traits, and decay rates, then builds a scenario, and runs the marsh equilibrium model over that scenario tracking individual sediment mass cohorts including minderal, live root, slow decay organic mater pool and a fast organic matter pool.
#' @param startYear an integer, year in form YYYY, the start year of the scenario 
#' @param endYear an integer, year in form YYYY, the end year of the scenario  
#' @param relSeaLevelRiseInit a numeric, initial rate of relative sea-level rise
#' @param rslrTotal a numeric, total relative sea-level rise over the course of the scanario
#' @param initElv a numeric, the initial elevation of the marsh at the start of the scenario
#' @param meanSeaLevel a numeric or a vector, Mean Sea Level at the start of the scenario, or a vector of Mean Sea Levels the same length as the number of years in a scenario
#' @param meanSeaLevelDatum a numeric, Mean Sea level over the last datum period
#' @param meanHighWater a numeric, Mean High Water level over the last datum period
#' @param meanHighHighWater a numeric (optional), Mean Higher High Water level over the last datum period
#' @param meanHighHighWaterSpring a numeric (optional), Mean Higher High Spring Tide Water level over the last datum period
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
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
#' @return a list of data frames, including the annualized summaries, and mapped cohorts tracked for every year of the simulation.
#' @export
runMemWithCohorts <- function(startYear, endYear=startYear+99, relSeaLevelRiseInit, rslrTotal, initElv,
                              meanSeaLevel, meanSeaLevelDatum=meanSeaLevel[1], meanHighWater, meanHighHighWater=NA, meanHighHighWaterSpring=NA, suspendedSediment, lunarNodalAmp,
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
  
  # Make sure tidyverse is there
  ##TODO move this to the package description
  require(tidyverse, quietly = TRUE)
  
  # Build scenario curve
  scenario <- buildScenarioCurve(startYear=startYear, endYear=endYear, meanSeaLevel=meanSeaLevel, 
                                 relSeaLevelRiseInit=relSeaLevelRiseInit, rslrTotal=rslrTotal, suspendedSediment=suspendedSediment)
  
  # add high tides
  scenario <- buildHighTideScenario(scenario, meanSeaLevelDatum=meanSeaLevelDatum, 
                                    meanHighWaterDatum=meanHighWater, meanHighHighWaterDatum=meanHighHighWater, meanHighHighWaterSpringDatum=meanHighHighWaterSpring, 
                                    lunarNodalAmp = lunarNodalAmp)
  
  # Add blank colums for attributes we will add later
  scenario$surfaceElevation <- as.numeric(rep(NA, nrow(scenario)))
  scenario$aboveground_biomass <- as.numeric(rep(NA, nrow(scenario)))
  scenario$belowground_biomass <- as.numeric(rep(NA, nrow(scenario)))
  scenario$mineral <- as.numeric(rep(NA, nrow(scenario)))

  # Set initial conditions
  # Calculate initial z star
    #TODO update the names to conform to new standards
  initialConditions <- determineInitialCohorts(initElv=initElv,
                                               MSL=MSL, MHW=MHW, MHHW=MHHW, MHHWS=MHHWS, ssc=ssc,
                                               bMax=bMax, zVegMin=zVegMin, zVegMax=zVegMax, zVegPeak=zVegPeak, 
                                               plantElevationType=plantElevationType,
                                               rootToShoot=rootToShoot, rootTurnover=rootTurnover, rootDepthMax=rootDepthMax, shape=shape,
                                               omDecayRate=omDecayRate, recalcitrantFrac=recalcitrantFrac, settlingVelocity=settlingVelocity,
                                               omPackingDensity=omPackingDensity, mineralPackingDensity=mineralPackingDensity,
                                               rootPackingDensity=omPackingDensity)
  cohorts <- initialConditions[[1]]

  
  # Add initial conditions to annual time step tracker
  scenario$surfaceElevation[1] <- initElv
  scenario$aboveground_biomass[1] <- initialConditions[[2]]
  scenario$belowground_biomass[1] <- initialConditions[[3]]
  scenario$mineral[1] <- initialConditions[[4]]
  
  cohorts$year <- rep(scenario$years[1], nrow(cohorts))
  
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
  
  # Fourth, add one cohort for each year in the scenario
  # Iterate through scenario table
  for (i in 2:nrow(scenario)) {
    
    # Calculate surface elevation relative to datum
    surfaceElvZStar <- zToZstar(z=scenario$surfaceElevation[i-1], meanHighWater=scenario$meanHighWater[i], meanSeaLevel=scenario$meanSeaLevel[i])
    
    # Calculate dynamic above ground
    dynamicAgb <- predictedBiomass(z=surfaceElvZStar, bMax = bMax, zVegMax = zStarVegMax, zVegMin = zStarVegMin, zVegPeak = zStarVegPeak)
    
    # Initial Below Ground Biomass
    dynamicBgb <- dynamicAgb * rootToShoot
    
    # Calculate Mineral pool
    dynamicMineralPool <- deliverSedimentFlexibly(z=scenario$surfaceElevation[i-1], suspendedSediment=scenario$suspendedSediment[i], 
                                                  meanSeaLevel=scenario$meanSeaLevel[i], meanHighWater=scenario$meanHighWater[i], meanHighHighWater = scenario$meanHighHighWater[i], 
                                                  meanHighHighWaterSpring = scenario$meanHighHighWaterSpring[i], settlingVelocity=settlingVelocity)
    
   
    # Calculate dynamic sediment deliver

    cohorts <- addCohort(cohorts, totalRootMassPerArea=dynamicBgb, rootDepthMax=rootDepthMax, 
                         rootTurnover = rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                         rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                         packing=list(organic=omPackingDensity, mineral=mineralPackingDensity), 
                         rootDensity=rootPackingDensity, shape=shape,
                         mineralInput = dynamicMineralPool)
    
    cohorts$year <- rep(scenario$years[i], nrow(cohorts))
    
    scenario$surfaceElevation[i] <- profileBottomElv + max(cohorts$layer_bottom, na.rm=T)
    scenario$aboveground_biomass[i] <- dynamicAgb
    scenario$belowground_biomass[i] <- dynamicBgb
    scenario$mineral[i] <- dynamicMineralPool

    # add initial set of cohorts
    cohortsNewRowMax <- cohortsNewRowMin + nrow(cohorts) - 1
    trackCohorts[cohortsNewRowMin:cohortsNewRowMax,] <- cohorts
    
    # create variables to keep track of cohorts added to the full cohort tracker
    cohortsNewRowMin <- cohortsNewRowMin + nrow(cohorts) + 1
  }
  
  # Remove NA values from cohorts
  trackCohorts <- trimCohorts(trackCohorts)
  
  # Return annual time steps and full set of cohorts
  outputsList <- list(annualTimeSteps = scenario, cohorts = trackCohorts)
  
  # If a core year is specificed return a core too.
  if (! is.na(coreYear)) {
    # Filter only to cohorts in the core year
    cohortsInCoreYear <- trackCohorts %>% 
      dplyr::filter(year==coreYear) %>%
      dplyr::select(-year)
    
    # Convert cohorts to age-depth profile
    coreYearAgeDepth <- convertProfileAgeToDepth(ageCohort=cohortsInCoreYear,
                                                  layerTop=coreMins,
                                                  layerBottom=coreMaxs)
    coreYearAgeDepth <- coreYearAgeDepth %>%  dplyr::filter(complete.cases(.)) %>%
      dplyr::mutate(om_fraction = (fast_OM+slow_OM+root_mass)/(fast_OM+slow_OM+root_mass+mineral),
             dry_bulk_density = (1/(
               (om_fraction/omPackingDensity) +
                 ((1-om_fraction)/mineralPackingDensity)
             )
             )
      )
    
    # Add to the list of outputs
    outputsList$core <- coreYearAgeDepth
  }

  return(outputsList)
}