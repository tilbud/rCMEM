#' Build high tides into scenario table
#'
#' This function takes a sea-level scenario table and elevations for 1 to 3 tidal cycles and outputs a scenario table with annualized high tidel levels.  
#' @param scenarioCurve a data frame, including year, mean sea level, and suspended sediment concentration, with each row being a year in the scenario
#' @param MSL0 a numeric, Mean Sea Level over a tidal datum period
#' @param MHW0 a numeric, Mean High Water level over the last datum period
#' @param MHHW0 a numeric, Mean Higher High Water level over the last datum period
#' @param MHHWS0 a numeric, Mean Higher High Spring Tide Water level over the last datum period
#' @param lunarNodalAmp the amplitude of the 18-year lunar nodal cycle
#' 
#' @return a data frame, including the sea-level and suspended sediment concentraiton scenario inputted, with annual high tide datum(s) added
#' @export
buildHighTideScenario <- function(scenarioCurve, MSL0=scenarioCurve$MSL[1], MHW0, MHHW0, MHHWS0, lunarNodalAmp) {
  
  # In create a MHW and add it to the scenario
  scenarioCurve$MHW <- predictLunarNodalCycle(year = scenarioCurve$years, MSL= scenarioCurve$MSL, 
                                             MSL0 = MSL0, floodElv=MHW0, 
                                             lunarNodalAmp=lunarNodalAmp)
  
  # If MHHW and MHHWS are arguements add them to the scenario table too
  if (!missing(MHHW0) & !missing(MHHWS0)) {
    scenarioCurve$MHHW <- predictLunarNodalCycle(year = scenarioCurve$years, MSL=scenarioCurve$MSL,
                                                MSL0 = MSL0, floodElv=MHHW0, 
                                                lunarNodalAmp=lunarNodalAmp)
    scenarioCurve$MHHWS <- predictLunarNodalCycle(year = scenarioCurve$years, MSL= scenarioCurve$MSL, 
                                                 MSL0 = MSL0, floodElv=MHHWS0, 
                                                 lunarNodalAmp=lunarNodalAmp)
  }
  return(scenarioCurve)
}