#' Add another cohort to the profile
#'
#' @param massPools 
#' @param consts 
#' @param dt_yr 
#' @param ssc 
#' @param meanTidalHeight 
#' @param rootTurnover 
#' @param rootOmFrac 
#' @param omDecayRate 
#' @param totalRootMass_per_area 
#' @param rootDepthMax 
#'
#' @return
#' @export
#'
#' @examples
addCohort <- function(massPools=data.frame(age=0, fast_OM=0, slow_OM=0, mineral=0, root_mass=0,
                                           layer_top=0, layer_bottom=0),
                      ssc, meanTidalHeight, #mineral inputs
                      rootTurnover, rootOmFrac, omDecayRate, #decay paraemters
                      totalRootMass_per_area, rootDepthMax, #layer thickness parameters
                      consts, dt_yr=1 ){
  if(!all(c('age', 'fast_OM', 'slow_OM', 'mineral', 'layer_top', 'layer_bottom') %in% names(massPools))){
    stop('Badly named massPools')
  }
  
  if(!all(c('packing') %in% names(consts))){
    stop('Can not find expected constants.')
  }
  
  packing <- consts$packing
  
  if(!all(c('organic', 'mineral') %in% names(packing))){
    stop('Can not find expected packing densities.')
  }
  
  if(!all(c('fast', 'slow') %in% names(rootOmFrac))){
    stop('Can not find expected root fraction splits.')
  }
  
  if(!all(c('fast', 'slow') %in% names(omDecayRate))){
    stop('Can not find expected organic matter decay rates.')
  }
  
  #decay the OM pools and add dead roots
  ans <- massPools %>% 
    dplyr::select(age, fast_OM, slow_OM, mineral, root_mass, layer_top, layer_bottom) %>%
    dplyr::mutate(age = age + dt_yr) %>%
    dplyr::mutate(fast_OM = fast_OM + 
             root_mass * rootOmFrac$fast * rootTurnover * dt_yr -
             fast_OM * omDecayRate$fast * dt_yr,
           slow_OM = slow_OM + 
             root_mass * rootOmFrac$slow * rootTurnover * dt_yr -
             slow_OM * omDecayRate$slow * dt_yr) %>%
    dplyr::select(-root_mass) #remove root_mass since it needs to be recacluated
  
  
  ans <- dplyr::bind_rows(data.frame(age = 0, fast_OM = 0, slow_OM = 0, 
                              mineral = sedimentInputs(ssc = ssc, 
                                                meanTidalHeight = meanTidalHeight, 
                                                consts=consts) * dt_yr),
                   ans) %>% #add sediments to the top
    dplyr::arrange(age) %>% #make sure things are sorted by age
    #calculate cumulative volumne of each pool
    dplyr::mutate(cumCohortVol = cumsum( (fast_OM + slow_OM)*packing$organic + mineral*packing$mineral )) %>%
  #calculate depth profile
    dplyr::mutate(layer_bottom = depthOfNotRootVolume(nonRootVolume = cumCohortVol,
                                                      totalRootMass_per_area =totalRootMass_per_area, 
                                                      rootDepthMax = rootDepthMax,
                                                      consts = consts)) %>%
    dplyr::mutate(layer_top = c(0, layer_bottom[-length(layer_bottom)])) %>%
    dplyr::mutate(root_mass = massLiveRoots(layerBottom = layer_bottom,
                                            layerTop = layer_top,
                                            totalRootMass_per_area =totalRootMass_per_area, 
                                            rootDepthMax = rootDepthMax,
                                            consts = consts))
  
  return(ans)
}