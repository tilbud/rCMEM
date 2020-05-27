#' Build scenario sea-level curve
#'
#' This function builds an annualied set of MEM inputs including sea-level rise and suspended sediment concentrations. 
#' @param startYear an integer, year in form YYYY, the start year of the scenario 
#' @param endYear an integer, year in form YYYY, the end year of the scenario  
#' @param MSL a numeric or a vector of numbers, either indicating mean sea-level at the start of scenario, or mean sea-level at each year of the sceanrio
#' @param rslr0 a numeric, initial rate of relative sea-level rise
#' @param rslrTotal a numeric, total relative sea-level rise over the course of the scanario
#' @param ssc a numeric or a vector of numbers, either average annual suspended sediment concentration, or a vector of annual suspended sediment concentration for each year of the scenario
#'
#' @references Sweet, W. V., Kopp, R. E., Weaver, C. P., Obeysekera, J., Horton, R. M., Thieler, E. R., & Zervas, C. (2017). Global and regional sea level rise scenarios for the United States.
#' @references Intergovernmental Panel on Climate Change (2013). Summary for policymakers, in Climate Change 2013: The Physical Science Basis, edited by T. F. Stocker, D. Qin, G.-K. Plattner, M. Tignor, S. K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex, and P. Midgley, pp. 3–29, Cambridge Univ. Press, Cambridge, U. K.
#' 
#' @return a data frame including columns for year, sea-level, and suspended sediment concentration, and rows for each year in the scenario
#' @export
buildScenarioCurve <- function(startYear, endYear=startYear+99, MSL, 
                               rslr0=0.3, rslrTotal=100, ssc) {
  
  # Create a sequence of the number of years
  years <- startYear:endYear
  nYearsOfSim <- length(years) # Put this in args. Replace this with yearEnd.
  
  scenario <- data.frame(index = 1:nYearsOfSim, years = years,
                         MSL = rep(NA, nYearsOfSim),
                         ssc = rep(NA, nYearsOfSim))
  
  # Build the Mean Sea Level Scenario
  # If the input only specifies an initial Mean Sea Level at time = 0 ... 
  if (length(MSL) == 1) {
    # ... create a sea-level rise scenario based on total sea-level rise and inital relative sea-level rise rate
    # Probably could use a stop here in case there are invalid entries of rslr0 and slrTotal
    scenario$MSL[1] <- MSL
    
    # IPCC, 2013 and Sweet et al., 2017
    # MSL(t) = MSL(0) + At + Bt^2 - where t is years from time zero
    # A=[MSL(1)-MSL(0)] – B  
    # Where [MSL(1)-MSL(0)] is the SLR at the start as determined from NOAA trend analyses 
    # Coefficient B is the acceleration term 
    # B = {[MSL(T)-MSL(0)] / T - [MSL(1)-MSL0] } / (T-1) where T - length of the simulation
    
    B <- ((rslrTotal)/nYearsOfSim - rslr0) / (nYearsOfSim-1)
    A <- rslr0 - B
    
    scenario$MSL[2:nYearsOfSim] <- scenario$MSL[1] + A*scenario$index[2:nYearsOfSim] + B*scenario$index[2:nYearsOfSim]^2
    
  } else if (length(MSL) == length(years)) {
    
    # If the user enters in a vector of MSL that is equal to the number of years in the simulation.
    scenario$MSL <- MSL
    
  } else {
    stop("RSLR input is incorrect. Either enter a value at for the beginning and ending of the scenario, or a vector of RSLR one for each year of the scenario.")
  }
  
  # Add suspended sediment concentration
  if (length(ssc) == 1) {
    # If ssc is a single value.
    scenario$ssc <- rep(ssc, nYearsOfSim)
  } else if (length(ssc) == nYearsOfSim) {
    # If ssc is a vector in equal length to the number of years of the simulation.
    scenario$ssc <- ssc
  } else {
    stop("SSC input is incorrect. Either enter a single average value, or a vector equal in length to the number of years in the scenario.")
  }
  return(scenario)
  # This is a start. We could also allow flexibility to change biomass and sediment addition functions later on to incorporate things like conversions of emergent to forested wetland.
}