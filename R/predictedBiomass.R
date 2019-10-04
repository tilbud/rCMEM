#' Predict aboveground biomass
#'
#'  Function takes elevation (z) as an input and outputs biomass.
#' Maximum growing biomass (bmax), upper and lower elevation limits
#' (zVegMin and zVegMax) are mandatory parameters
#' The elevation of peak biomass is an optional parameter.
#' This is a version of the function that can create skewed parabolic distributions.
#'
#' @param z marsh elevation
#' @param bMax maximum biomass
#' @param zVegMax upper elevation of biomass limit
#' @param zVegMin lower elevation of biomass limit
#' @param zVegPeak (optional) elevation of peak biomass
#'
#' @return a numerical value for the aboveground biomass
#' @export
#'
#' @examples
#' predictedBiomass(0, 2500, 3, -1)
predictedBiomass <- function(z=0, bMax=2500, zVegMax=3, zVegMin=-1, zVegPeak=NA) {

  
  # Stop the function if there are invalid parameters
  if ( (bMax < 0) | # negative peak biomss
       # or elevations that don't make sense
       (zVegMax <= max(zVegMin, zVegPeak, na.rm=T)) |
       ((zVegMin >= min(zVegMax, zVegPeak, na.rm=T)))
       ) {
    stop("invalid biomass parameters")
  }
  
  # If there is no peak elevation for vegetation parabola is symmetric.
  if (is.na(zVegPeak)) { 
    
    # Elevation of the veg. peak is halfway between the min and max limits
    zVegPeak<-(zVegMax+zVegMin)/2
  
    # From bmax, min, and max elevation limits, solve for parameters of a parabola.
    a <- -((-zVegMin * bMax - zVegMax * bMax) / ((zVegMin - zVegPeak) * (-zVegMax + zVegPeak)))
    b <- -(bMax / ((zVegMin - zVegPeak) * (-zVegMax + zVegPeak)))
    c <- (zVegMin * zVegMax * bMax) / ((zVegMin - zVegPeak) * (zVegMax - zVegPeak))
    
    # Apply parabolic function to elevation to calulate above ground biomass.
    agb <- a*z + b*z^2 + c
  
  # If elevation of peak biomass is specified the curve can be more flexible ...
  } else {
    # ... but we need to split it into two curves.
    
    # For the curve applied on the 'upper end' create a new minimum elevation
    # mirroring  upper elevation limit accross the peak biomass elevation.
    zVegMin_up <- zVegPeak-((zVegMax-zVegPeak)) 
    
    # For the curve applied at the 'lower end', same. Create new maximum elevation 
    # tolerance mirroring, lower elevation limit accros the peak biomass elevation.
    zVegMax_low <-zVegPeak+((zVegPeak-zVegMin))
    
    # Solve for the parameters of the upper curve.
    a_up <- -((-zVegMin_up * bMax - zVegMax * bMax) / ((zVegMin_up - zVegPeak) * (-zVegMax + zVegPeak)))
    b_up <- -(bMax / ((zVegMin_up - zVegPeak) * (-zVegMax + zVegPeak)))
    c_up <- (zVegMin_up * zVegMax * bMax) / ((zVegMin_up - zVegPeak) * (zVegMax - zVegPeak))
    
    # Solve for the parametrs of the lower curve.
    a_low <- -((-zVegMin * bMax - zVegMax_low * bMax) / ((zVegMin - zVegPeak) * (-zVegMax_low + zVegPeak)))
    b_low <- -(bMax / ((zVegMin - zVegPeak) * (-zVegMax_low + zVegPeak)))
    c_low <- (zVegMin * zVegMax_low * bMax) / ((zVegMin - zVegPeak) * (zVegMax_low - zVegPeak))
    
    # If elevation is above the specified peak biomass elevation, apply the upper curve,
    # if it's under apply the lower curve
    agb <- ifelse(z>zVegPeak, a_up*z + b_up*z^2 + c_up, a_low*z + b_low*z^2 + c_low)
    
  }
  
  # Recast any negative values as 0, since biomass can't be negative.
  agb <- ifelse(agb>0,agb,0) 
  return(agb) # Return biomass as a function of elevation.
  
}