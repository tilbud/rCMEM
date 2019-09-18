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
#' @param soilLength unit length of soil area, generally 1
#' @param soilWidth unit width of soil area, generally 1.
#' @param shape flag of root shape, only \code{linear} is implimented.
#' @param ... 
#' @param rootDepthMax an integer that is the maximum root depth, generally in cm
#'
#' @return an array of the associated live root mass for each layer specified.
#' @export
#'
#' @examples
#' massLiveRoots(layerBottom = 1:10, layerTop = 0:9, 
#'               totalRootMass_per_area = 0.3, rootDepthMax=5,
#'               soilLength=1, soilWidth=1, shape='linear')
massLiveRoots <- function(layerBottom, layerTop, 
                          totalRootMass_per_area, 
                          rootDepthMax,
                          soilLength=1, soilWidth=1, 
                          shape,
                          expDecayRate_perMaxDepth = log(0.05),
                          ...){
  
  totalRootMass <- soilLength*soilWidth*totalRootMass_per_area
  
  if (totalRootMass == 0) {
    rootMass <- rep(0, length(layerBottom))
    return(rootMass)
  } else {
    
    layerBottom[is.na(layerBottom)] <- 0
    layerTop[is.na(layerTop)] <- 0
    
    if(any(layerBottom < layerTop)){
      #print(data.frame(layerBottom, layerTop))
      stop('Bad layer definition.')
    }
    
    ##reset the layers that are beyond the root inputs to have 0 depths
    layerBottom[layerBottom > rootDepthMax] <- rootDepthMax
    layerTop[layerTop > rootDepthMax] <- rootDepthMax
    
    if(shape == 'linear'){
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
      
    }else{
      if(shape == 'exponential'){
        # root mass with depth x: r(x) = a * exp(b * x) - m
        # Let R(x) be the antiderivative: R(x) = a / b * exp(b * x) - m * x + C0
        # Assume
        # r( x = rootDepthMax) => m = a * exp(b * rootDepthMax)
        # Total root biomass is integral of r from 0 to rootDepthMax
        #                 TotalMass  = a / b * exp(b*rootDepthMax) - m * rootDepthMax - a / b
        #                            = a / b * exp(b*rootDepthMax) -  a * exp(b * rootDepthMax) * rootDepthMax - a / b
        #  a = TotalMass * [1 / b * exp(b * rootDepthMax) -  exp(b * rootDepthMax) * rootDepthMax - 1 / b]^-1
        #
        # mass between two layers is the integral of r from x1 to x2
        #         mass = (a / b * exp(b * x2) - m * x2) - (a / b * exp(b * x1) - m * x1)
        #         mass = a / b * (exp(b * x2)-exp(b * x1)) + m *(x1 - x2)
        #layerTop <- 0; layerBottom = 1; expDecayRate_perRootDepthMax <- log(0.05); totalRootMass <- 0.2; rootDepthMax <- 30
        b <- expDecayRate_perMaxDepth / rootDepthMax #convert from total profile to per cm
        a <- totalRootMass * (1 / b * exp(b * rootDepthMax) -  
                                exp(b * rootDepthMax) * rootDepthMax - 1 / b)^-1
        m <- a * exp(b * rootDepthMax)
        rootMass <- a / b * (exp(b * layerBottom)-exp(b * layerTop)) +  m*(layerTop - layerBottom)
        
      }else{
        
        stop(paste('Unknown shape specified:', shape))
      }
    }
    return(rootMass)
  }
  
  
  
}
