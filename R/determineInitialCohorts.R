#' Determine the initial conditions for a set of cohorts
#' @param initElv a numeric, the initial elevation of the marsh at the start of the scenario
#' @param meanSeaLevel a numeric, Mean Sea Level
#' @param meanHighWater a numeric, Mean High Water level
#' @param meanHighHighWater a numeric (optional), Mean Higher High Water level
#' @param meanHighHighWaterSpring a numeric (optional), Mean Higher High Spring Tide Water level
#' @param suspendedSediment a numeric, suspended sediment concentration of the water column
#' @param captureRate a numeric, the number of times a water column will clear per tidal cycle
#' @param nFloods a numeric, the number of tidal flooding events per year
#' @param floodTime.fn a function, specify the method used to calculate flooding time per tidal cycle
#' @param bMax a numeric or vector of numerics, maximum biomass
#' @param zVegMax a numeric or vector of numerics, upper elevation of biomass limit
#' @param zVegMin a numeric or vector of numerics, lower elevation of biomass limit
#' @param zVegPeak a numeric or vector of numerics, (optional) elevation of peak biomass
#' @param plantElevationType character, either "orthometric" or "dimensionless", specifying elevation reference of the vegetation growing elevations
#' @param rootToShoot a numeric or vector of numerics, root to shoot ratio
#' @param rootTurnover a numeric or vector of numerics, below ground biomass annual turnover rate
#' @param abovegroundTurnover (optional) a numeric or vector of numerics, aboveground biomass annual turnover rate
#' @param speciesCode (optional) a character or vector of characters, species names or codes associated with biological inputs
#' @param rootDepthMax a numeric or vector of numerics, maximum (95\%) rooting depth
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
#' 
#' @return a data frame of mass pools representing the initial conditions for a simulation
#' @export
determineInitialCohorts <- function(initElv,
                                 meanSeaLevel, meanHighWater, meanHighHighWater=NA, 
                                 meanHighHighWaterSpring=NA, suspendedSediment,
                                 nFloods = 705.79, floodTime.fn = floodTimeLinear,
                                 bMax, zVegMin, zVegMax, zVegPeak, plantElevationType,
                                 rootToShoot, rootTurnover, rootDepthMax, shape="linear",
                                 abovegroundTurnover=NA, speciesCode=NA,
                                 omDecayRate, recalcitrantFrac, captureRate,
                                 omPackingDensity=0.085, mineralPackingDensity=1.99,
                                 rootPackingDensity=omPackingDensity,
                                 initialCohorts=NA,
                                 uplandCohorts=NA,
                                 supertidalCohorts=NA,
                                 supertidalSedimentInput=NA,
                                 ...) {
  
  # If initialCohorts is defined as an input then it overrides all other arguements.
  if (is.data.frame(initialCohorts)) {
    # If it does,
    # Return initial cohorts
    cohorts <- initialCohorts
    bio_table <- data.frame(speciesCode=NA,
                            rootToShoot=NA,
                            rootTurnover=NA,
                            abovegroundTurnover=NA,
                            rootDepthMax=NA,
                            aboveground_biomass=NA,
                            belowground_biomass=NA)
    initSediment <- NA
  } else { # If initial cohorts are not supplied
    
    # Convert real growing elevations to dimensionless growing elevations
    if (! plantElevationType %in% c("dimensionless", "zStar", "Z*", "zstar")) {
      zStarVegMin <- convertZToZstar(zVegMin, meanHighWater, meanSeaLevel)
      zStarVegMax <- convertZToZstar(zVegMax, meanHighWater, meanSeaLevel)
      zStarVegPeak <- convertZToZstar(zVegPeak, meanHighWater, meanSeaLevel)
    } else {
      zStarVegMin <- zVegMin
      zStarVegMax <- zVegMax
      zStarVegPeak <- zVegPeak
    }
    
    # Convert dimensionless plant growing elevations to real growing elevations
    if (plantElevationType %in% c("dimensionless", "zStar", "Z*", "zstar")) {
      zVegMin <- convertZStarToZ(zVegMin, meanHighWater, meanSeaLevel)
      zVegMax <- convertZStarToZ(zVegMax, meanHighWater, meanSeaLevel)
      zVegPeak <- convertZStarToZ(zVegPeak, meanHighWater, meanSeaLevel)
    }
    
    # Set initial conditions
    # Calculate initial z star
    initElvStar <- convertZToZstar(z=initElv, meanSeaLevel=meanSeaLevel, meanHighWater=meanHighWater)
    
    # Initial Above Ground Biomass
    bio_table <- runMultiSpeciesBiomass(initElvStar, bMax = bMax, zVegMax = zStarVegMax, 
                                        zVegMin = zStarVegMin, zVegPeak = zStarVegPeak,
                                        rootToShoot=rootToShoot,
                                        rootTurnover=rootTurnover, 
                                        abovegroundTurnover=abovegroundTurnover, 
                                        rootDepthMax=rootDepthMax, speciesCode=speciesCode)
    
    # If elevation is lower than highest tide provided, and lower than maximum growing elevation
    # Generate 1 m or more of sediment given equilibrium conditions
    if ((initElv <= max(meanHighWater, meanHighHighWater, meanHighHighWaterSpring, na.rm=T)) & (initElv <= max(zVegMax))) {

      # Initial Sediment
      initSediment <- deliverSediment(z=initElv, suspendedSediment=suspendedSediment, nFloods=nFloods,
                                      meanSeaLevel=meanSeaLevel, meanHighWater=meanHighWater,
                                      meanHighHighWater=meanHighHighWater, meanHighHighWaterSpring=meanHighHighWaterSpring,
                                      captureRate=captureRate,
                                      floodTime.fn=floodTime.fn)
      
      # Run initial conditions to equilibrium
      cohorts <- runToEquilibrium(totalRootMassPerArea=bio_table$belowground_biomass[1], rootDepthMax=bio_table$rootDepthMax[1],
                                  rootTurnover=bio_table$rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                                  rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                                  packing=list(organic=omPackingDensity, mineral=mineralPackingDensity),
                                  rootDensity=rootPackingDensity, shape=shape, 
                                  mineralInput = initSediment,
                                  minDepth = round(max(rootDepthMax)+0.5))
      
    } else if ((initElv >= max(meanHighWater, meanHighHighWater, meanHighHighWaterSpring, na.rm=T)) & (initElv <= max(zVegMax))) { 
      # If elevation is greater than highest tide provided, but lower than maximum growing elevation
      # Then form super-tidal peat
      
      # Check to see if supertidal peat is defined as an input
      if (is.data.frame(supertidalCohorts)) {
        # If it is, than pass the input staight to the output
        cohorts <- supertidalCohorts
        bio_table <- data.frame(speciesCode=NA,
                                rootToShoot=NA,
                                rootTurnover=NA,
                                abovegroundTurnover=NA,
                                rootDepthMax=NA,
                                aboveground_biomass=NA,
                                belowground_biomass=NA)
        initSediment <- NA
      } else {
        # If supertidalSedimentInput is defined
        if (is.data.frame(supertidalSedimentInput)) {
          # Run initial conditions to equilibrium
          initSediment <- supertidalSedimentInput
          cohorts <- runToEquilibrium(totalRootMassPerArea=bio_table$belowground_biomass[1], rootDepthMax=bio_table$rootDepthMax[1],
                                      rootTurnover=bio_table$rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                                      rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                                      packing=list(organic=omPackingDensity, mineral=mineralPackingDensity),
                                      rootDensity=rootPackingDensity, shape=shape, 
                                      mineralInput = initSediment,
                                      minDepth = round(max(rootDepthMax)+0.5))
        } else {
          # If not, come up with a set with a column of of peat generated with biomass inputs,
          # and any assumed 0 sediment input.
          
          # Initial Sediment, an arbitrary low number
          # Here I use the annual sediment delivered 1 cm below the highest tide line
          initSediment <- deliverSediment(z=max(meanHighWater, meanHighHighWater, meanHighHighWaterSpring, na.rm=T)-1, 
                                          suspendedSediment=suspendedSediment, 
                                          meanSeaLevel=meanSeaLevel, 
                                          meanHighWater=meanHighWater, 
                                          meanHighHighWater=meanHighHighWater, 
                                          meanHighHighWaterSpring=meanHighHighWaterSpring,
                                          nFloods = nFloods,
                                          captureRate=captureRate,
                                          floodTime.fn = floodTime.fn)
          
          cohorts <- runToEquilibrium(totalRootMassPerArea=bio_table$belowground_biomass[1], rootDepthMax=bio_table$rootDepthMax[1],
                                      rootTurnover=bio_table$rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                                      rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                                      packing=list(organic=omPackingDensity, mineral=mineralPackingDensity),
                                      rootDensity=rootPackingDensity, shape=shape, 
                                      mineralInput = initSediment,
                                      minDepth = round(max(rootDepthMax)+0.5))
        }
      }
    } else if ((initElv >= max(meanHighWater, meanHighHighWater, meanHighHighWaterSpring, na.rm=T)) & (initElv > zVegMax)) {
      
      # If elevation is greater than maximum growing elevation
      # Then assign it an upland soil
      # If an upland soil is provided use it.
      if (! is.na(uplandCohorts)) {
        cohorts <- uplandCohorts
        bio_table <- data.frame(speciesCode=NA,
                                rootToShoot=NA,
                                rootTurnover=NA,
                                abovegroundTurnover=NA,
                                rootDepthMax=NA,
                                aboveground_biomass=NA,
                                belowground_biomass=NA)
        initSediment <- NA
      } else {
        # If not assign it an arbitrary 50% organic matter soil
        cohorts <- data.frame(age=rep(0, round(rootDepthMax+0.6)),
                              fast_OM=rep(0, round(rootDepthMax+0.6)),
                              slow_OM=rep(0.5*(1/(0.5/mineralPackingDensity+0.5/omPackingDensity)), round(rootDepthMax+0.6)),
                              respired_OM=rep(0, round(rootDepthMax+0.6)),
                              mineral=rep(0.5*(1/(0.5/mineralPackingDensity+0.5/omPackingDensity)),round(rootDepthMax+0.6)),
                              root_mass=rep(0,round(rootDepthMax+0.6)),
                              layer_top=0:((round(rootDepthMax+0.6)-1)),
                              layer_bottom=1:round(rootDepthMax+0.6)) %>% 
          dplyr::mutate(cumCohortVol = cumsum(layer_bottom-layer_top))
        
        bio_table <- data.frame(speciesCode=NA,
                                rootToShoot=NA,
                                rootTurnover=NA,
                                abovegroundTurnover=NA,
                                rootDepthMax=NA,
                                aboveground_biomass=NA,
                                belowground_biomass=NA)
        initSediment <- NA
      }
    } else {
      stop("Elevations are invalid for creating initial cohorts.")
    }
  }
  
  # Check to make sure it has the right column names,
  # If it does then return it,
  # If not throw an error message
  if (! all(c("age", "fast_OM", "slow_OM", 
              "mineral", "root_mass", 
              "layer_top", "layer_bottom", "cumCohortVol") %in% names(cohorts))) {
    missing <- paste(c("age", "fast_OM", "slow_OM", 
                       "mineral", "root_mass", 
                       "layer_top", "layer_bottom")[!c("age", "fast_OM", "slow_OM", 
                                                       "mineral", "root_mass", 
                                                       "layer_top", "layer_bottom", "cumCohortVol")%in%names(cohorts)],
                     collapse = ", ")
    stop(paste("Initial cohorts table is missing ", missing, ".", sep=""))
  }
  
  return(list(cohorts,
              bio_table,
              initSediment))
}
