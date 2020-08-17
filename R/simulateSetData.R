#' Simulate SET Data
#' 
#' This will take MEM outputs and simulate sediment elevation table and marker horizon data.
#' 
#' @param cohorts a dataframe, one output from runMemWithCohorts 
#' @param scenario a dataframe, a second output from runMemWithCohorts
#' @param markerHorizonYear an integer, the year a marker horizon was lain down, that were using to simulate an accretion measurement 
#' 
#' @return a dataframe of net elevation changes and accretion rates relative to year
#' 
#' @export
simulateSetData <- function(cohorts, scenario, markerHorizonYear=NA) {
  # For the scenario, 
  # For each set of cohorts measure the minimum depth of the cohort nearest to the horizon age without going older
  # Accretion rate is mimium depth / (time t - marker horizon year)
  scenario <- scenario %>% 
    dplyr::mutate(netElevationChange = surfaceElevation - lag(surfaceElevation))
  
  if (! is.na(markerHorizonYear)) {
    cohortsMh <- cohorts %>% 
      dplyr::filter(complete.cases(.)) %>% 
      dplyr::group_by(year) %>% 
      dplyr::filter(year >= markerHorizonYear) %>%
      dplyr::mutate(index = length(age):1)
    
    mhYear <- cohortsMh$index[1]
    
    cohortsMh <- cohortsMh %>% 
      dplyr::filter(index == mhYear) %>%
      dplyr::mutate(accretionRate = layer_top / age) %>% 
      dplyr::select(year, accretionRate) # %>% 
      # dplyr::rename(years = year)
    
    cohortsMh$accretionRate[1] <- 0 
    
    scenario <- scenario %>% 
      dplyr::left_join(cohortsMh, by="year")
  }
  return(scenario)
}