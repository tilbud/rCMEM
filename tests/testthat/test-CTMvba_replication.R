context("Repeat vba CTM results")

test_that("run3 test",{
  ##Comment out when not in dev
  load('data/CTM_vba_testcase.RData')
  purrr::walk(list.files('R', full.names = TRUE), source)
  
  parm <- setNames(as.list(testCases$parameters.df$Run3), testCases$parameters.df$Sheet)
  
  profile_byAge <- data.frame(age=0, fast_OM=0, slow_OM=0, mineral=0, root_mass=0,
                              layer_top=0, layer_bottom=0)

for(ii in 1:210){
  profile_byAge <- addCohort(massPools = profile_byAge,
                               rootTurnover = parm$BGTR, 
                               rootOmFrac = list(fast=(1-parm$kr), slow=parm$kr),
                               omDecayRate = list(fast=parm$OMDR, slow=0),
                               packing = list(root = parm$k1, #density of roots g-C/cm3
                                     organic = parm$k1,  # Organic Self Packing Densities: k1 g-C/cm3
                                     mineral = parm$k2),
                               
                               mineralInput_g_per_yr.fn = sedimentInputs,
                               ssc = parm$SSC,##Play with a correction factor
                               meanTidalHeight = parm$D, 
                               
                               massLiveRoots.fn = massLiveRoots,
                               totalRootMass_per_area = parm$RT/1e4,
                               rootDepthMax=parm$Dm,
                               
                               depthOfNotRootVolume.fn = depthOfNotRootVolume,
                               rootDensity = parm$k1)
}
  endingProfile_byAge <- profile_byAge %>% 
    filter(is.finite(age), (layer_bottom-layer_top) > 0, age > 0) %>%
    mutate(layer_bottom = layer_bottom - layer_top[1]) %>%
    mutate(layer_top = layer_top -layer_top[1]) 
  
  endingProfile_byDepth <- convertProfile_AgeToDepth(endingProfile_byAge, 
                             layerTop = c(0, testCases$Run3$`Core Section or Depth (cm)`[-nrow(testCases$Run3)]),
                             layerBottom = testCases$Run3$`Core Section or Depth (cm)`)
  
  compareReady <- endingProfile_byDepth %>%
    mutate("Core Section or Depth (cm)" = layer_bottom,
           "years of input" = input_yrs,
           "Calc'd volume (cm3)" = layer_bottom-layer_top,
           "Horizon Age (yr)" = age,
           ##conver from g to to g/m2 to 
           "Root Biomass (g/m2)" = (root_mass)*1e4,
           "Labile OM (g/m2)" = fast_OM/((layer_bottom-layer_top)*1e-2*1e-2),
           "Refrac Biomass (g/m2)" = slow_OM/((layer_bottom-layer_top)*1e-2*1e-2),
           "Mineral (g/cm2)" = mineral/((layer_bottom-layer_top)*1e-2*1e-2)) %>%
    select(`Core Section or Depth (cm)`:`Mineral (g/cm2)`) %>%
    left_join(testCases$Run3 %>% select(`Core Section or Depth (cm)`:`Mineral (g/cm2)`),
              by='Core Section or Depth (cm)', suffix = c('___R', '___VBA')) %>%
    gather(key='variable_model', value='value', -`Core Section or Depth (cm)`)%>%
    tidyr::separate(col=variable_model, into=c('varable', 'model'), sep='___') 
  
  ggplot(compareReady, aes(x=`Core Section or Depth (cm)`, y=value, color=model)) +
    geom_line() +
    facet_wrap(~varable, scales='free')
  
  
})