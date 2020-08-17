#' Simulate a sediment core from a set of cohorts
#' @param cohorts a data frame
#' @param coreYear an integer, year in form YYYY, specify a year to simulate taking a sediment core
#' @param coreDepth an integer, depth, specify a depth to simulate coring to and assume 1 cm sampling intervals
#' @param coreMins a vector of sampling depth minimums to simulate coring subsables, this is an alternative to depth, and 1cm increments
#' @param coreMaxs a vector of sampling depth maximums to simulate coring subsables, this is an alternative to depth, and 1cm increments
#' @param omToOcParams a list of parameters defining a linear or quadratic relationship between fraction organic matter and fraction carbon
#' @param omPackingDensity a numeric, the bulk density of pure organic matter
#' @param mineralPackingDensity a numeric, the bulk density of pure mineral matter
#' @param rootPackingDensity a numeric, the bulk density of pure root matter
#' 
#' @return a dataframe with variables simulated from a soil core
#' @export
simulateSoilCore <- function(cohorts, coreYear, coreDepth=100, 
                             coreMaxs=1:coreDepth, coreMins=coreMaxs-1,
                             omToOcParams = list(B0=0, B1=0.48),
                             omPackingDensity=0.085, mineralPackingDensity=1.99,
                             rootPackingDensity=omPackingDensity) {
  
  # Filter only to cohorts in the core year
  cohortsInCoreYear <- cohorts %>% 
    dplyr::filter(year==coreYear) %>%
    dplyr::select(-year)
  
  # If the profile is shorter than the intended core ...
  if (max(coreMaxs) > max(cohortsInCoreYear$layer_bottom, na.rm=T)) {
    coreMins <- coreMins[coreMins< max(cohortsInCoreYear$layer_bottom, na.rm=T)]
    coreMaxs <- coreMaxs[1:length(coreMins)]
    coreMaxs[length(coreMaxs)] <- max(cohortsInCoreYear$layer_bottom, na.rm=T)
  }
  # 
  
  # To convert OM to OC
  # If parmeter list is 2 long then simple linear correlation
  if (length(omToOcParams) == 2) {
    omToOc <- function(om, B0=omToOcParams$B0, B1=omToOcParams$B1) {return(B0 + om*B1)}
  } else if (length(omToOcParams) == 3) {
    # If parameter list is 3 long, then it's quadratic
    omToOc <- function(om, B0=omToOcParams$B, B1=omToOcParams$B1,
                       B2=omToOcParams$B2) {return(B0 + om*B1 + om^2*B2)}
  } else {
    # If something else then trip an error message
    stop("Invalid number of organic matter to organic carbon conversion parameters,")
  }
  
  # Convert cohorts to age-depth profile
  coreYearAgeDepth <- convertProfileAgeToDepth(ageCohort=cohortsInCoreYear,
                                                layerTop=coreMins,
                                                layerBottom=coreMaxs)
  
  coreYearAgeDepth <- coreYearAgeDepth %>%  
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(om_fraction = (fast_OM+slow_OM+root_mass)/(fast_OM+slow_OM+root_mass+mineral),
                  dry_bulk_density = (1/(
                    (om_fraction/omPackingDensity) +
                      ((1-om_fraction)/mineralPackingDensity))
                    ),
                  oc_fraction = omToOc(om_fraction)
                  )
  
  # Add to the list of outputs
  return(coreYearAgeDepth) 
}