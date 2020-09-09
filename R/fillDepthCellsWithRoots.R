#' Calculates a root mass depth profile
#' 
#'
#' It will take any root biomass and maximum rooting depth, 
#' as well as vectors of minimum and maximum depths.
#' We use predefined equations if rootShape is set to linear or exponential
#' It can also take a custom function of x, as long as you define parameters in
#' a custom parameters list.
#' The custom function listed here is a gaussian root depth distribution.
#'
#' @param depthMins a vector of numerics, minimum of depth intervals that the funciton will fill with roots
#' @param depthMaxs a vector of numerics, maximum of depth intervals that the funciton will fill with roots
#' @param totalRootBmass a numeric, total root biomass
#' @param rootShape a character, either linear, exponential, or custom, defining the shape of the root profile
#' @param rootDepthMax a numeric, maximum depth of root zone
#' @param customFunction a function (optional), if root shape is custom, then use a customized function 
#' @param customParams a list (optional), any other parameters required by a custom root shape function
#'
#' @return a data frame with minimum and maximum depths and and root biomass within those depths intervals
#' @export
#'
fillDepthCellsWithRoots <- function(depthMins = 0:149, 
                                    depthMaxs = 1:150,
                                    totalRootBmass = 3000,
                                    rootShape = 'linear',
                                    rootDepthMax = 30,
                                    customFunction = function(x) { 
                                      exp(-((x-maxRootDepth)^2/lambdaRoot^2)) 
                                    },
                                    customParams = list(lambdaRoot = 33,
                                                        maxRootDepth = 30)
) {
  

  
  # Calculate root mass for each min to max depth inteval.
  if (rootShape == "linear") { # if we're calculating a linear root mass shape
    
    # Calculate intercept of linear function
    rootMassInt <- 2 * totalRootBmass / rootDepthMax
    
    # If the depth inteval is above the rooting depth
    # calculate root mass as integral between min to max depth.
    # If it's below rootin depth, root mass is 0.
    rootMassInCell <- ifelse(depthMins < rootDepthMax,
                             rootMassInt * (depthMaxs - depthMins) -
                               rootMassInt * (depthMaxs ^ 2 - depthMins ^ 2) /
                               (2 * rootDepthMax),
                             0)
    
  } else if (rootShape == "exponential") { 
    
    # if we're calculating an exponential root mass shape
    # Define an asymtote that is a small number 
    rootBmassAsymtote <- log(0.05) / rootDepthMax
    
    # Calculate intercept of the exponential function
    rootMassInt <- -0.95 * totalRootBmass * rootBmassAsymtote / 
      (1 - exp(rootBmassAsymtote * rootDepthMax))
    
    #  calculate root mass as integral between min to max in two steps.
    e1 <- exp(rootBmassAsymtote * depthMaxs) # step 1
    
    # step 2
    rootMassInCell <- (rootMassInt / rootBmassAsymtote) * 
      (e1 - exp(depthMins * rootBmassAsymtote))
    
  } else if (rootShape == "custom") { # If another custom function,
    
    # Go through the list and make sure each parameter is an object in memory
    for (i in 1:length(customParams)) {
      paramName <- names(customParams[i])
      assign(paramName, customParams[[paramName]])
    }
    
    # Iterate through depth series and integrate function at min and max depth
    rootDensity <- c() # blank vector for storing outputs
    for (i in 1:length(depthMaxs)) {
      rootMassCell_i <-integrate(customFunction, depthMins[i], depthMaxs[i]
      )$value
      rootDensity <- c(rootDensity, rootMassCell_i)
    }
    
    # normalize and mutliply by total biomass
    rootMassInCell <- rootDensity / sum(rootDensity) * totalRootBmass
    
  }
  
  # Return a data frame with depth intervals and roots in cells
  return(data.frame(depthMin = depthMins, 
                    depthMax = depthMaxs, 
                    rootBiomass_gPerM2 = rootMassInCell))
}