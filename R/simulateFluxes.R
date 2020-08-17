#' Simulate Greenhouse Gas Fluxes
#' 
#' This function takes outputs from runMemWithCohorts and uses them to simulate several annualized greenhouse gas flux variables 
#' 
#' @param cohorts
#' @param scenario
#' @param omToOcParams
#' @param salinity
#' @param salThreshold
#' @param salRate
#' @param maxCH4effect
#' 
#' @return list of 2 data frames, one the input cohorts table with depthwise fluxes added, the second the annualized scenario table with annualized fluxes added
#' @export
simulateFluxes <- function(cohorts, scenario,
                           omToOcParams = list(B0=0, B1=0.48),
                           salinity = 32,
                           salThreshold = 18,
                           salRate = 0.5,
                           maxCh4ToCo2 = 1/5) {
  
  # To convert OM to OC
  # If parmeter list is 2 long then simple linear correlation
  if (length(omToOcParams) == 2) {
    omToOc <- function(om, B0=omToOcParams$B0, B1=omToOcParams$B1) {return(B0 + om*B1)}
  } else if (length(omToOcParams) == 3) {
    # If parameter list is 3 long, then it's quadratic
    omToOc <- function(om, B0=omToOcParams$B0, B1=omToOcParams$B1,
                       B2=omToOcParams$B2) {return(B0 + om*B1 + om^2*B2)}
  } else {
    # If something else then trip an error message
    stop("Invalid number of organic matter to organic carbon conversion parameters,")
  }
  
  # Constant used throughout
  hoursInTidalCycle <- 12.42
  
  # If one single salinity value is entered add it to the annual time steps table all
  # If multiple salinity values are entered, add to the annual time step table,
  # as long as they are the same lenght. 
  # If mutliple values are entered, but not the same number of years of the scenario, trip an error.
  if (length(salinity == 1) | length(salinity) == nrow(scenario)) {
    scenario <- dplyr::mutate(scenario, salinity = salinity)
  } else {
    stop("Invalid entry for salinities. Use either one single value, or one for each year of the simulation.")
  }
  
  # Create a table of tidal cycles
  if ((! any(is.na(scenario$meanHighWater))) & any(is.na(scenario$meanHighHighWater)) & any(is.na(scenario$meanHighHighWaterSpring))) {
    # If only MHW is included
    highTidesPerYear <- 352.657 + 352.657
    tidalCycles <- data.frame(datum = c("meanHighWater"),
                              nTides = c(highTidesPerYear),
                              stringsAsFactors = F)
    
  } else if ((! any(is.na(scenario$meanHighWater))) & (! any(is.na(scenario$meanHighHighWater))) & any(is.na(scenario$meanHighHighWaterSpring))) {
    # If MHW and MHHW are included
    highTidesPerYear <- 352.657
    higherHighTidesPerYear <- 352.657
    
    tidalCycles <- data.frame(datum = c("meanHighWater", "meanHighHighWater"),
                              nTides = c(highTidesPerYear, higherHighTidesPerYear),
                              stringsAsFactors = F)
    
  } else if (all(!is.na(c(scenario$meanHighWater, scenario$meanHighHighWater, scenario$meanHighHighWaterSpring)))) {
    # If MHHWS, MHHW, and MHW all included
    highTidesPerYear <- 352.657
    higherHighTidesPerYear <- 352.657 - 24.72
    springTidesPerYear <- 24.72
    
    tidalCycles <- data.frame(datum = c("meanHighWater", "meanHighHighWater", "meanHighHighWaterSpring"),
                              nTides = c(highTidesPerYear, higherHighTidesPerYear, springTidesPerYear),
                              stringsAsFactors = F)
  } else {
    # Else stop
    stop("Invalid tidal datums")
  }
  
  timeStepsRelevantData <- scenario %>% 
    dplyr::select(year, meanSeaLevel, matches("meanHighWater|meanHighHighWater|meanHighHighWaterSpring"), surfaceElevation, salinity) %>% 
    # TODO rename years year in upstream code so we don't have to do it here
    # rename(year=years) %>% 
    tidyr::gather(key = "datum", value = "datumHigh",
                  -year, -meanSeaLevel, -surfaceElevation, -salinity) %>% 
    dplyr::mutate(datumLow = meanSeaLevel-(datumHigh-meanSeaLevel)) %>% 
    dplyr::left_join(tidalCycles, by="datum") %>% 
    dplyr::arrange(year, datum) # not really necessary
  
  C_to_CH4 <- 16.04 / 12.01
  C_to_CO2 <- 44.01 / 12.01
  
  ghgFluxesCohorts <- cohorts %>% 
    dplyr::select(year, age, respired_OM, layer_top, layer_bottom) %>% 
    # Join the Surface Elevation, MSL, MHW, MHHW, and MHHWS, 
    # and salinity data from annual time steps to cohorts
    dplyr::full_join(timeStepsRelevantData, by="year") %>% 
    # Calculate fractional inundation time for each cohort
    dplyr::mutate(layer_mid = layer_bottom - ((layer_bottom-layer_top)/2),
                  z = surfaceElevation - layer_mid,
                  floodTime = ifelse(z > datumHigh,
                                     0, floodTimeFromDatum(z = z, datumHigh = datumHigh,
                                                           datumLow = datumLow) * nTides)) %>%
    dplyr::select(year, age, layer_top, layer_bottom, layer_mid, respired_OM, salinity, floodTime) %>% 
    dplyr::group_by(year, age, layer_top, layer_bottom, layer_mid, respired_OM, salinity) %>% 
    dplyr::summarise(floodTime = sum(floodTime)) %>% 
    dplyr::mutate(floodFraction = floodTime / 8760,
                  respired_C = omToOc(respired_OM),
                  salinity_effect = maxCh4ToCo2 / (1+exp(salRate*(salinity-salThreshold))),
                  respired_C_CH4 = respired_C * floodFraction * salinity_effect,
                  respired_C_CO2 = respired_C - respired_C_CH4,
                  respired_CH4 = respired_C_CH4 * C_to_CH4,
                  respired_CO2 = respired_C_CO2 * C_to_CO2)
  
  # Each cohorts outgoing C in CH4 is a function of
  # Inundation time and whether sulfates turn methane production on or off
  # Each cohort's outgoing C in CO2 is the inverse
  # Convert C to CH4
  # Convert C to CO2
  
  ghgFluxesAnnual <- ghgFluxesCohorts %>% 
    dplyr::group_by(year) %>% 
    # Summarise each years CH4 emission
    # Summarise each year's CO2 respiration
    dplyr::summarise(respired_CH4 = sum(respired_CH4),
                     respired_CO2 = sum(respired_CO2))
  
  co2_removal <- cohorts %>% 
    dplyr::group_by(year) %>% 
    dplyr::summarise(slow_OM = sum(slow_OM)) %>% 
    dplyr::mutate(slow_OM_change = slow_OM - lag(slow_OM),
                  sequestered_C = -omToOc(slow_OM_change),
                  sequestered_CO2 = sequestered_C * C_to_CO2) %>% 
    select(year, sequestered_CO2)
  
  ghgFluxesScenario <- scenario %>% 
    dplyr::mutate(belowgroundNPP = omToOc(belowground_biomass * rootTurnover) * C_to_CO2,
                  abovegroundNPP = omToOc(aboveground_biomass * abovegroundTurnover) * C_to_CO2,
                  totalNPP = abovegroundNPP + belowgroundNPP) %>% 
    left_join(ghgFluxesAnnual, by=c("year"))
  
  return(list(ghgFluxesCohorts, ghgFluxesScenario))
  
}
