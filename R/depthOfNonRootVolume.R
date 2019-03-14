#' Depth of the given soil volume
#' 
#' Given a shape and mass for the roots, calculate the depth of a specific non-root volume.
#'
#' @param nonRootVolume a numerical value that is the volume of interest, generally in cm3.
#' @param totalRootMass_per_area an integer that is the total mass per area of the roots, generally in g cm-3
#' @param rootDepthMax an integer that is the maximum root depth, generally in cm
#' @param parms placeholder currently NULL
#' @param consts a list of constants that must include numerics representing the soilLength (cm) and soilWidth (cm) for the area of interest, a flag shape, and a list called packing that includes a value of the root packing (or density). 
#'
#' @return
#' @export
#'
#' @examples
depthOfNotRootVolume <- function(nonRootVolume, 
                                 totalRootMass_per_area, 
                                 rootDepthMax,
                                 parms=NULL, consts){
 
  if(!all(c('soilLength', 'soilWidth', 'shape', 'packing') %in% names(consts))){
    stop('Can not find all constants.')
  }
  soilLength <- consts$soilLength
  soilWidth <- consts$soilWidth
  rootDensity <- consts$packing$root
  shape <- consts$shape
  
  ####
  totalRootMass <- soilLength*soilWidth*totalRootMass_per_area
  totalRootVolume <- totalRootMass/rootDensity
  
  if(totalRootVolume > soilLength*soilWidth*rootDepthMax){
    stop('Bad root volume')
  }
  
  if(shape == 'linear'){
    rootWidth <- totalRootVolume*2/(rootDepthMax*soilLength)

    #nonRootVolume = ((rootWidth/rootDepthMax*depth^2)/2 + depth*(soilWidth-rootWidth))*soilLength
    #            0 = rootWidth / (2*rootDepthMax) * depth ^2 + 
    #                   (soilWidth-rootWidth) * depth - nonRootVolume/soilLength ##solve for depth
    coef1 <- rootWidth / (2*rootDepthMax)
    coef2 <- soilWidth-rootWidth
    coef3 <- -nonRootVolume/soilLength
    ansDepth <- (-coef2 + sqrt(coef2^2-4*coef1*coef3))/(2*coef1)
    
    #correct for beyond root zone
    behondRootZone <- nonRootVolume > soilLength*soilWidth*rootDepthMax - totalRootVolume
    ansDepth[behondRootZone] <- (rootDepthMax +
                            (nonRootVolume - soilLength*soilWidth*rootDepthMax + totalRootVolume ) /
                              (soilLength*soilWidth)) [behondRootZone] #treat as a square
    return(ansDepth)
  }else{
    stop('Unknown shape specified')
  }
}
