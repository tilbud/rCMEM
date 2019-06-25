context("Repeat vba CTM results")

test_that("run3 test",{
  ##Comment out when not in dev
  load('data/CTM_vba_testcase.RData')
  purrr::walk(list.files('R', full.names = TRUE), source)
  
  parm <- setNames(as.list(testCases$parameters.df$Run3), testCases$parameters.df$Sheet)
  
  startingProfile <- data.frame(age=0, fast_OM=0, slow_OM=0, mineral=0, root_mass=0,
                              layer_top=0, layer_bottom=0)

for(ii in 1:183){
  startingProfile <- addCohort(massPools = startingProfile,
                               ssc = parm$SSC, 
                               meanTidalHeight = parm$D, 
                               rootTurnover = parm$BGTR, 
                               rootOmFrac = (1-parm$kr),
                               omDecayRate = defaultParms$omDecayRate,
                               totalRootMass_per_area = defaultParms$totalRootBiomass, 
                               rootDepthMax = defaultParms$rootDepthMax,
                               consts = defaultConsts)
}
})