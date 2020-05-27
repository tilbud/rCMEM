#' Run Marsh Equilibrium Model with Cohorts
#'
#' This function takes an initial elevation, average annual suspended sediment concentrtion and tidal properties, plant traits, and decay rates, then builds a scenario, and runs the marsh equilibrium model over that scenario tracking individual sediment mass cohorts including minderal, live root, slow decay organic mater pool and a fast organic matter pool.
#' @param startYear an integer, year in form YYYY, the start year of the scenario 
#' @param endYear an integer, year in form YYYY, the end year of the scenario  
#' @param rslrT1 a numeric, initial rate of relative sea-level rise
#' @param rslrTotal a numeric, total relative sea-level rise over the course of the scanario
#' @param initElv a numeric, the initial elevation of the marsh at the start of the scenario
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
#' @param coreYear an integer, year in form YYYY, (optional) specify a year to simulate taking a sediment core
#' @param coreDepth an integer, depth, (optional) specify a depth to simulate coring to and assume 1 cm sampling intervals
#' @param coreMins a vector of sampling depth minimums (optional) to simulate coring subsables, this is an alternative to depth, and 1cm increments
#' @param coreMaxs a vector of sampling depth maximums (optional) to simulate coring subsables, this is an alternative to depth, and 1cm increments
#' @param ...
#' 
#' @return a list of data frames, including the annualized summaries, mapped cohorts tracked for every years, and if core year is specified, a core.
#' @export
runMemWithCohorts <- function(startYear, endYear=startYear+99, rslrT1, rslrTotal, initElv,
                              MSL, MSL0=MSL[1], MHW, MHHW=NA, MHHWS=NA, ssc, lunarNodalAmp,
                              bMax, zVegMin, zVegMax, zVegPeak, plantElevationType,
                              rootToShoot, rootTurnover, rootDepthMax, shape="linear",
                              omDecayRate, recalcitrantFrac, settlingVelocity,
                              omPackingDensity=0.085, mineralPackingDensity=1.99,
                              rootPackingDensity=omPackingDensity,
                              coreYear=NA, coreDepth=100, coreMaxs=1:coreDepth, coreMins=coreMaxs-1,
                              ...) {
  
  # Make sure tidyverse is there
  require(tidyverse, quietly = TRUE)
  
  # Build scenario curve
  scenario <- buildScenarioCurve(startYear=startYear, endYear=endYear, MSL=MSL, 
                                 rslr0=rslrT1, rslrTotal=rslrTotal, ssc=ssc)
  
  # add high tides
  scenario <- buildHighTideScenario(scenario, MSL0=MSL0, 
                                    MHW0=MHW, MHHW0=MHHW, MHHWS0=MHHWS, 
                                    lunarNodalAmp = lunarNodalAmp)
  
  # Add blank colums for attributes we will add later
  scenario$surfaceElevation <- as.numeric(rep(NA, nrow(scenario)))
  scenario$biomass <- as.numeric(rep(NA, nrow(scenario)))
  scenario$mineral <- as.numeric(rep(NA, nrow(scenario)))
  
  # Convert dimensionless plant growing elevations to real growing elevations
  if (! plantElevationType %in% c("dimensionless", "zStar", "Z*", "zstar")) {
    zStarVegMin <- zToZstar(zVegMin, MHW, MSL[1])
    zStarVegMax <- zToZstar(zVegMax, MHW, MSL[1])
    zStarVegPeak <- zToZstar(zVegPeak, MHW, MSL[1])
  } else {
    zStarVegMin <- zVegMin
    zStarVegMax <- zVegMax
    zStarVegPeak <- zVegPeak
  }
  
  # Set initial conditions
  # Calculate initial z star
  initElvStar <- zToZstar(z=initElv, MSL=scenario$MSL[1], MHW=scenario$MHW[1])
  
  # Initial Above Ground Biomass
  initAgb <- predictedBiomass(z=initElvStar, bMax = bMax, zVegMax = zStarVegMax, 
                              zVegMin = zStarVegMin, zVegPeak = zStarVegPeak)
  
  # Initial Below Ground Biomass
  initBgb <- initAgb * rootToShoot
  
  # Initial Sediment
  initSediment <- deliverSedimentFlexibly(z=initElv, ssc=scenario$ssc[1], 
                                          MSL=scenario$MSL[1], MHW=scenario$MHW[1], 
                                          MHHW = scenario$MHHW[1], MHHWS = scenario$MHHWS[1],
                                          settlingVelocity=settlingVelocity)
  
  # Run initial conditions to equilibrium
  cohorts <- runToEquilibrium(totalRootMass_per_area=initBgb, rootDepthMax=rootDepthMax,
                              rootTurnover=rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                              rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                              packing=list(organic=omPackingDensity, mineral=mineralPackingDensity),
                              rootDensity=rootPackingDensity, shape=shape, 
                              mineralInput_g_per_yr = initSediment)
  
  # Add initial conditions to annual time step tracker
  scenario$surfaceElevation[1] <- initElv
  scenario$mineral[1] <- initSediment
  scenario$biomass[1] <- initBgb
  
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
    surfaceElvZStar <- zToZstar(z=scenario$surfaceElevation[i-1], MHW=scenario$MHW[i], MSL=scenario$MSL[i])
    
    # Calculate dynamic above ground
    dynamicAgb <- predictedBiomass(z=surfaceElvZStar, bMax = bMax, zVegMax = zStarVegMax, zVegMin = zStarVegMin, zVegPeak = zStarVegPeak)
    
    # Initial Below Ground Biomass
    dynamicBgb <- dynamicAgb * rootToShoot
    
    # Calculate Mineral pool
    dynamicMineralPool <- deliverSedimentFlexibly(z=scenario$surfaceElevation[i-1], ssc=scenario$ssc[i], 
                                                  MSL=scenario$MSL[i], MHW=scenario$MHW[i], MHHW = scenario$MHHW[i], 
                                                  MHHWS = scenario$MHHWS[i], settlingVelocity=settlingVelocity)
    
   
    # Calculate dynamic sediment deliver
    cohorts <- addCohort(cohorts, totalRootMass_per_area=dynamicBgb, rootDepthMax=rootDepthMax, 
                         rootTurnover = rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                         rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                         packing=list(organic=omPackingDensity, mineral=mineralPackingDensity), 
                         rootDensity=rootPackingDensity, shape=shape,
                         mineralInput_g_per_yr = dynamicMineralPool)
    
    cohorts$year <- rep(scenario$years[i], nrow(cohorts))
    
    scenario$surfaceElevation[i] <- profileBottomElv + max(cohorts$layer_bottom, na.rm=T)
    scenario$biomass[i] <- dynamicBgb
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
    coreYearAgeDepth <- convertProfile_AgeToDepth(ageCohort=cohortsInCoreYear,
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