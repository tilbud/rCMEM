#' Run to equilbrium with proscribed depth to mean high tide
#' 
#' 
#'
#' @param parms a list of parameters that are passed to addCohort individually
#' @param consts a list of constants that are passed to addCohort as a list
#' @param minAge run the model for atleast these many years
#' @param maxAge do not run the model for longer then these many years
#' @param recordEvolution a boolean flag to record how the equalibrium profile evolves
#' @param relTol stop the evolution if the relative difference between the two oldest cohorts is less then this
#' @param absTol stop the evolution if the absolute difference between the two oldest cohorts is less then this
#'
#' @return a data frame by age cohort with the layer top and bottom depths, age, and mass of the fast organic, slow organic, mineral, and root biomass pools.
#' @export
#'
#' @examples
runToEquilibrium <- function(parms, consts, 
                             minAge = 50, maxAge = 200,
                             recordEvolution = FALSE,
                             relTol = 1e-6, absTol = 1e-8){
  
  #initalize things to empty
  cohortProfile <- data.frame(age=0, fast_OM=0, slow_OM=0, mineral=0, root_mass=0,
                              layer_top=NA, layer_bottom=NA)
  
  if(recordEvolution){
    record.ls <- list(cohortProfile)
  }
  for(ii in 2:maxAge){
    if(recordEvolution){
      record.ls[[sprintf('Yr%d', ii)]] <- cohortProfile
    }
    #oldCohort <- cohortProfile
    cohortProfile <- addCohort(massPools = cohortProfile,
                               ssc = parms$ssc, 
                               meanTidalHeight = parms$depthBelowMHW, 
                               rootTurnover = parms$rootTurnover, 
                               rootOmFrac = parms$rootOmFrac,
                               omDecayRate = parms$omDecayRate,
                               totalRootMass_per_area = parms$totalRootBiomass, 
                               rootDepthMax = parms$rootDepthMax,
                               consts = consts)
    
    if(ii == 2){
      cohortProfile <- cohortProfile[1,] #chop off the 0 initalizaiton
    }
    ##have the last layer OM pools stabilized?
    if(ii > minAge){
      if((abs(diff(cohortProfile$fast_OM[ii-c(1:2)] + cohortProfile$slow_OM[ii-c(1:2)])) < absTol |
         abs(diff(cohortProfile$fast_OM[ii-c(1:2)] + cohortProfile$slow_OM[ii-c(1:2)] ) /
             (cohortProfile$fast_OM[ii-1] + cohortProfile$slow_OM[ii-1])) < relTol)){
        break
      }
    }
  }
  
  if(recordEvolution){
    record.ls
  }else{
    return(cohortProfile)
  }
}