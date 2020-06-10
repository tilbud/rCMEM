#' Determine the initial conditions for a set of cohorts
#' @param initElv a numeric, the initial elevation of the marsh at the start of the scenario
#' @param MSL a numeric, Mean Sea Level over the last datum period
#' @param MHW a numeric, Mean High Water level over the last datum period
#' @param MHHW a numeric (optional), Mean Higher High Water level over the last datum period
#' @param MHHWS a numeric (optional), Mean Higher High Spring Tide Water level over the last datum period
#' @param ssc a numeric, suspended sediment concentration of the water column
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
#' 
#' @return a data frame of mass pools representing the initial conditions for a simulation
#' @export
determineInitialCohorts <- function(initElv,
                                 MSL, MHW, MHHW=NA, MHHWS=NA, ssc,
                                 bMax, zVegMin, zVegMax, zVegPeak, plantElevationType,
                                 rootToShoot, rootTurnover, rootDepthMax, shape="linear",
                                 omDecayRate, recalcitrantFrac, settlingVelocity,
                                 omPackingDensity=0.085, mineralPackingDensity=1.99,
                                 rootPackingDensity=omPackingDensity,
                                 initialCohorts=NA,
                                 uplandCohorts=NA,
                                 supertidalCohorts=NA,
                                 supertidalSedimentInput=NA) {
  
  # If initialCohorts is defined as an input then it overrides all other arguements.
  if (! is.na(initialCohorts)) {
    # If it does,
    # Return initial cohorts
    cohorts <- initialCohorts
    initBgb <- NA
    initSediment <- NA
  } else { # If initial cohorts are not supplied
    
    # Convert real growing elevations to dimensionless growing elevations
    if (! plantElevationType %in% c("dimensionless", "zStar", "Z*", "zstar")) {
      zStarVegMin <- zToZstar(zVegMin, MHW, MSL)
      zStarVegMax <- zToZstar(zVegMax, MHW, MSL)
      zStarVegPeak <- zToZstar(zVegPeak, MHW, MSL)
    } else {
      zStarVegMin <- zVegMin
      zStarVegMax <- zVegMax
      zStarVegPeak <- zVegPeak
    }
    
    # Convert dimensionless plant growing elevations to real growing elevations
    if (plantElevationType %in% c("dimensionless", "zStar", "Z*", "zstar")) {
      zVegMin <- zStarToZ(zVegMin, MHW, MSL)
      zVegMax <- zStarToZ(zVegMax, MHW, MSL)
      zVegPeak <- zStarToZ(zVegPeak, MHW, MSL)
    }
    
    # Set initial conditions
    # Calculate initial z star
    initElvStar <- zToZstar(z=initElv, MSL=MSL, MHW=MHW)
    
    # Initial Above Ground Biomass
    initAgb <- predictedBiomass(z=initElvStar, bMax = bMax, zVegMax = zStarVegMax, 
                                zVegMin = zStarVegMin, zVegPeak = zStarVegPeak)
    
    # Initial Below Ground Biomass
    initBgb <- initAgb * rootToShoot
    
    # If elevation is lower than highest tide provided, and lower than maximum growing elevation
    # Generate 1 m or more of sediment given equilibrium conditions
    if ((initElv <= max(MHW, MHHW, MHHWS, na.rm=T)) & (initElv <= zVegMax)) {

      # Initial Sediment
      initSediment <- deliverSedimentFlexibly(z=initElv, ssc=ssc, 
                                              MSL=MSL, MHW=MHW, 
                                              MHHW=MHHW, MHHWS=MHHWS,
                                              settlingVelocity=settlingVelocity)
      
      # Run initial conditions to equilibrium
      cohorts <- runToEquilibrium(totalRootMass_per_area=initBgb, rootDepthMax=rootDepthMax,
                                  rootTurnover=rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                                  rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                                  packing=list(organic=omPackingDensity, mineral=mineralPackingDensity),
                                  rootDensity=rootPackingDensity, shape=shape, 
                                  mineralInput_g_per_yr = initSediment,
                                  minDepth = round(rootDepthMax+0.5))
      
    } else if ((initElv >= max(MHW, MHHW, MHHWS, na.rm=T)) & (initElv <= zVegMax)) { 
      # If elevation is greater than highest tide provided, but lower than maximum growing elevation
      # Then form super-tidal peat
      
      # Check to see if supertidal peat is defined as an input
      if (! is.na(initialCohorts)) {
        # If it is, than pass the input staight to the output
        cohorts <- supertidalCohorts
        initBgb <- NA
        initSediment <- NA
      } else {
        # If supertidalSedimentInput is defined
        if (! is.na(supertidalSedimentInput)) {
          # Run initial conditions to equilibrium
          initSediment <- supertidalSedimentInput
          cohorts <- runToEquilibrium(totalRootMass_per_area=initBgb, rootDepthMax=rootDepthMax,
                                      rootTurnover=rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                                      rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                                      packing=list(organic=omPackingDensity, mineral=mineralPackingDensity),
                                      rootDensity=rootPackingDensity, shape=shape, 
                                      mineralInput_g_per_yr = initSediment,
                                      minDepth = round(rootDepthMax+0.5))
        } else {
          # If not, come up with a set with a column of of peat generated with biomass inputs,
          # and any assumed 0 sediment input.
          
          # Initial Sediment, an arbitrary low number
          # Here I use 1/10 of the annual sediment delivered 1 cm below the highest tide line
          initSediment <- deliverSedimentFlexibly(z=max(MHW, MHHW, MHHWS, na.rm=T)-1, 
                                                  ssc=ssc, 
                                                  MSL=MSL, MHW=MHW, 
                                                  MHHW=MHHW, MHHWS=MHHWS,
                                                  settlingVelocity=settlingVelocity)
          
          cohorts <- runToEquilibrium(totalRootMass_per_area=initBgb, rootDepthMax=rootDepthMax,
                                      rootTurnover=rootTurnover, omDecayRate = list(fast=omDecayRate, slow=0),
                                      rootOmFrac=list(fast=1-recalcitrantFrac, slow=recalcitrantFrac),
                                      packing=list(organic=omPackingDensity, mineral=mineralPackingDensity),
                                      rootDensity=rootPackingDensity, shape=shape, 
                                      mineralInput_g_per_yr = initSediment,
                                      minDepth = round(rootDepthMax+0.5))
          cohorts <- cohorts %>% 
            dplyr::mutate(cumCohortVol = cumsum(layer_bottom-layer_top))
        }
      }
    } else if ((initElv >= max(MHW, MHHW, MHHWS, na.rm=T)) & (initElv > zVegMax)) {
      
      # If elevation is greater than maximum growing elevation
      # Then assign it an upland soil
      # If an upland soil is provided use it.
      if (! is.na(uplandCohorts)) {
        cohorts <- uplandCohorts
        initBgb <- NA
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
        
                              
        initBgb <- 0
        initSediment <- 0
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
  
  # Check to make sure the initial cohorts are deeper down than the maximum rooting depth
  #if (max(cohorts$layer_bottom) < rootDepthMax) {
  #  stop(paste("Initial cohort depth in less than maximum rooting depth. Either reenter initial cohorts, or increase the maxAge setting on runToEquilibrium."))
  #}
  
  return(list(cohorts,
              initAgb,
              initBgb,
              initSediment))
}
