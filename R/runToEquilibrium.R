#' Run to equilbrium with proscribed surface mineral input
#' 
#' Repeatedly add new cohorts to a soil profile until the oldest layers are not siginficatly different AND the oldest layer are above some minimum age, within some upper bound of soil age.
#'
#' @param minAge run the model for atleast these many years
#' @param maxAge do not run the model for longer then these many years
#' @param recordEvolution a boolean flag to record how the equalibrium profile evolves
#' @param relTol stop the evolution if the relative difference between the two oldest cohorts is less then this
#' @param absTol stop the evolution if the absolute difference between the two oldest cohorts is less then this
#' @param ... arguments to be passed to \code{addCohort}
#'
#' @return a data frame by age cohort with the layer top and bottom depths, age, and mass of the fast organic, slow organic, mineral, and root biomass pools OR a list that is a record of age cohorts
#' @export
#'
runToEquilibrium <- function(minAge = 50, maxAge = 200,
                             recordEvolution = FALSE,
                             relTol = 1e-6, absTol = 1e-8, ...){
  
  #initalize things to empty
  cohortProfile <- data.frame(age=NA, fast_OM=NA, slow_OM=NA, 
                              respired_OM=NA,
                              mineral=NA, root_mass=NA,
                              layer_top=NA, layer_bottom=NA)
  
  if(recordEvolution){
    record.ls <- list(cohortProfile)
  }
  for(ii in 2:maxAge){
    if(recordEvolution){
      record.ls[[sprintf('Yr%d', ii)]] <- cohortProfile
    }
    #oldCohort <- cohortProfile
    cohortProfile <- addCohort(massPools = cohortProfile, ...)
    
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
    return(record.ls)
  }else{
    return(cohortProfile)
  }
}