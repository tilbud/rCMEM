#' Mass of living roots
#'
#' This function returns the mass of the living roots between two layers. Given a mass per area of roots and a maximum rooting depth, this calculates the total mass between two layers with a specified distribution. Currently only a linear distribution is implimented. \cr
#' The linear algorithm is as follows:  \itemize{
#'    \item #mass_per_depth = slope * depth + intercept
#'    \item slope <- -2 * totalRootMass / (rootDepthMax^2)
#'    \item intercept <- 2 * totalRootMass / rootDepthMax
#'    \item #mass = integral(mass_per_depth, depth)
#'    \item rootMass <- intercept * (layerBottom-layerTop) + slope/2 * (layerBottom ^2-layerTop^2)
#'    }
#' 
#' @param layerBottom an array of depths from soil top to the bottom of the layer, generally in cm
#' @param layerTop an array of depths from soil top to the top of the layer, generally in cm
#' @param totalRootMass_per_area an integer that is the total mass per area of the roots, generally in g cm-3
#' @param rootDepthMax an integer that is the maximum root depth, generally in cm
#' @param parms currently a placeholder and not used
#' @param consts a list of constants including integers for the soilLength and soilWidth of the area of interest, as well as the flag for the shape of the root mass. See example for details.
#'
#' @return an array of the associated live root mass for each layer specified.
#' @export
#'
#' @examples
#' massLiveRoots(layerBottom = 1:10, layerTop = 0:9, 
#'               totalRootMass_per_area = 0.3, rootDepthMax=5,
#'               consts=list(soilLength=1, soilWidth=1, shape='linear'))
massLiveRoots <- function(layerBottom, layerTop, 
                          totalRootMass_per_area, 
                          rootDepthMax,
                          parms = NULL, consts){
  
  if(!all(c('soilLength', 'soilWidth', 'shape') %in% names(consts))){
    stop('Can not find all constants.')
  }
  
  totalRootMass <- consts$soilLength*consts$soilWidth*totalRootMass_per_area
  
  if (totalRootMass == 0) {
    rootMass <- rep(0, length(layerBottom))
    return(rootMass)
  } else {
    
    if(any(layerBottom < layerTop)){
    stop('Bad layer definition.')
  }
  
  ##reset the layers that are beyond the root inputs to have 0 depths
  layerBottom[layerBottom > rootDepthMax] <- rootDepthMax
  layerTop[layerTop > rootDepthMax] <- rootDepthMax
  
  if(consts$shape == 'linear'){
    #mass_per_depth(depth) = slope * depth + intercept; mass(depth) = slope/2*depth^2+intercept*depth
    #Given
    #   mass_per_depth(maxDepth) = 0  => intercept = -slope * maxDepth
    #   mass(maxDepth) = TotalMass => 
    #                                 TotalMass = slope/2*maxDepth^2 + intercept*maxDepth
    #                                 => slope = 2*TotalMass/maxDepth^2
    #                                    intercept = -2 * TotalMass/maxDepth
    #
    slope <- -2 * totalRootMass / (rootDepthMax^2)
    intercept <- 2 * totalRootMass / rootDepthMax
    #mass = integral(mass_per_depth, depth)
    rootMass <- intercept*(layerBottom-layerTop) + slope/2*(layerBottom ^2-layerTop^2)
    return(rootMass)
  }else{
 
    stop('Unknown shape specified')
  }
  }
    
  
  
}
