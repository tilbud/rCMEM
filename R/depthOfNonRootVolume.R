#' Depth of the given soil volume
#' 
#' Given a shape and mass for the roots, calculate the depth of a specific non-root volume.
#'
#' @param nonRootVolume a numerical value (possibly array) that is the volume of interest, generally in cm3.
#' @param totalRootMass_per_area an integer that is the total mass per area of the roots, generally in g cm-3
#' @param rootDepthMax an integer that is the maximum root depth, generally in cm
#' @param rootDensity a numeric that is the root density in g per cm3
#' @param soilLength a numeric of the unit length of soil volume, generally 1 cm
#' @param soilWidth a numeric of the unitl witdth of the soil colume, generally 1cm
#' @param shape a character of the shape of the root distribution, currently only 'linear' is implimented
#' @param ... 
#'
#' @return A numeric (possibly array) corresponding to the depth of the specificed non root volume
#' 
#' @export
depthOfNotRootVolume <- function(nonRootVolume, 
                                 totalRootMass_per_area, 
                                 rootDepthMax,
                                 rootDensity,
                                 soilLength = 1, soilWidth=1,
                                 shape='linear', ...){
  
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
