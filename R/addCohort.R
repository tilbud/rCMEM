#' Add another cohort to the soil profile
#'
#' This function adds a new cohort to a soil profile. Cohorts are constructed as follows:
#' \enumerate{
#'        \item Cohorts are aged by dt_yr incrament.
#'        \item Organic matter pools decay
#'        \item If there are non-zero mineral inputs, they are added to the top of the soil profile as a new cohort.
#'        \item The top and bottom of the cohort layer is recalculated.
#'        \item A new root profile is constructed based on the cohort layers.
#' }
#'
#' @param massPools A data frame containing the following colums: age, fast_OM, slow_OM, mineral, layer_top, layer_bottom, root_mass. Each row of this data frame represents a single cohort. To increase runtime in simulations, the profile is buffered with NA's at the top of the data frame, see code comments for details.
#' @param packing A list of numerics with mineral and organic arguments representing their respective packing densities. 
#' @param dt_yr A numeric of the time step in years. Typically set to 1 year.
#' @param rootTurnover A numeric as the root turnover time in faction per year.
#' @param rootOmFrac A list of numerics called fast and slow which is the allocation fraction between fast and slow organic matter pool respectively for dead roots.
#' @param omDecayRate A list of numerics called fast and slow which is the decay rate for the fast and slow organic matter pool in fraction per year.
#' @param massLiveRoots.fn A function that returns the mass of live roots for specified depth layers, must accept \code{layerBottom} and \code{layerTop} as arguments.
#' @param depthOfNonRootVolume.fn A function that returns the depth of a specified volumen of soil, must accept \code{nonRootVolume} as an argument.
#' @param mineralInput_g_per_yr.fn A function that returns the mineral input in grams per year (numeric)
#' @param ... arguments to be passed to the specified functions
#'
#' @return A new data frame similar to massPool with the added cohort and simulated change in the cohort pools.
#' @export
#' 
#'
addCohort <- function(massPools,
                      rootTurnover, rootOmFrac, omDecayRate, #decay paraemters
                      packing, #packing densities
                      mineralInput_g_per_yr.fn = sedimentInputs, 
                      mineralInput_g_per_yr = NA,
                      massLiveRoots.fn = massLiveRoots,
                      depthOfNonRootVolume.fn = depthOfNonRootVolume,
                      dt_yr=1, ...){
  
  #Sanity check the inputs
  if(!all(c('age', 'fast_OM', 'slow_OM', 'mineral', 'layer_top', 'layer_bottom', 'root_mass') %in% 
          names(massPools))){
    stop('Badly named massPools')
  }
  
  if(!all(c('organic', 'mineral') %in% names(packing))){
    stop('Can not find expected packing densities.')
  }
  
  if(!all(c('fast', 'slow') %in% names(rootOmFrac))){
    stop('Can not find expected root fraction splits.')
  }
  
  if(!all(c('fast', 'slow') %in% names(omDecayRate))){
    stop('Can not find expected organic matter decay rates.')
  }
  
  #copy over to an answer data frame
  ans <- massPools
  
  ans$age <- ans$age + dt_yr #age the cohorts
  
  # track respiration
  ans$respired_OM <- ans$fast_OM +
    (ans$root_mass * rootOmFrac$fast * rootTurnover * dt_yr) -
    ans$fast_OM * omDecayRate$fast * dt_yr
  
  #add and decay the organic matter
  ans$fast_OM <- ans$fast_OM + 
             ans$root_mass * rootOmFrac$fast * rootTurnover * dt_yr -
             ans$fast_OM * omDecayRate$fast * dt_yr
  
  ans$slow_OM <- ans$slow_OM + 
             ans$root_mass * rootOmFrac$slow * rootTurnover * dt_yr -
             ans$slow_OM * omDecayRate$slow * dt_yr
  

  
  # Check to see if mineral input is a static value or a function
  if (is.na(mineralInput_g_per_yr)) {
    mineralInput_g_per_yr <- mineralInput_g_per_yr.fn(...)
  }
  
  #if we are laying down a new cohort
  if(mineralInput_g_per_yr > 0){
    if(!any(is.na(ans$age))){ #if there aren't any empty cohort slots
      bufferAns <- ans
      bufferAns[TRUE] <- NA
      ans <- rbind(bufferAns, ans) #double the number of slots
    }
    
    #Take the last empty cohort slot
    ##Note that by buffering the profile with empty cohort slots we do not need to copy over the 
    ##...entire data frame when adding a single cohort. This increases runtime when this function is
    ##...embedided in interative runs.
    newCohortIndex <- max(which(is.na(ans$age)))
    
    #initalize the age and organic carbon pools to 0
    ans$age[newCohortIndex] <- 0
    ans$fast_OM[newCohortIndex] <- 0
    ans$slow_OM[newCohortIndex] <- 0
    ans$respired_OM[newCohortIndex] <- 0
    
    #lay down the new mineral inputs
    ans$mineral[newCohortIndex] <- mineralInput_g_per_yr * dt_yr
    
  }
  
  #calculate the volume of each cohort
  temp_Vol <- (ans$fast_OM + ans$slow_OM)/packing$organic + ans$mineral/packing$mineral
  temp_Vol[is.na(temp_Vol)] <- 0 #replace NA with 0 so we can calculate cumulative colume
  ans$cumCohortVol <- cumsum(temp_Vol)
  
  #calculate depth profile
  ans$layer_bottom <- depthOfNonRootVolume.fn(nonRootVolume = ans$cumCohortVol,
                                              massLiveRoots.fn=massLiveRoots.fn,
                                              soilLength=1, soilWidth=1,
                                              relTol = 1e-6,
                                              ...)
  ans$layer_top <- c(0, ans$layer_bottom[-length(ans$layer_bottom)])
  
  #recalculate the root mass
  ans$root_mass <- massLiveRoots.fn(layerBottom = ans$layer_bottom,
                                    layerTop = ans$layer_top, ...)
  
  return(ans)
}