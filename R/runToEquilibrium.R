runToEquilibrium <- function(cohortStep = addCohort, parms, consts, 
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
    cohortProfile <- cohortStep(massPools = cohortProfile,
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