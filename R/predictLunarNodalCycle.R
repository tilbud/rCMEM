#' Predict lunar nodal cycle
#'
#' This function builds a high tide level based on the 18.1 year lunar nodal cycle.
#' @param MSL0 a numeric, mean sea-level over the last complete 18 year datum period
#' @param MSL a vector, mean sea-level for each scenario we are predicting
#' @param floodElv a numeric, high tide level over the last complete 18 year datum period
#' @param lunarNodalAmp the amplitude of the 18 year lunar nodal cycle
#' 
#' @references Morris, J. T., Sundberg, K., & Hopkinson, C. S. (2013). Salt marsh primary production and its responses to relative sea level and nutrients in estuaries at Plum Island, Massachusetts, and North Inlet, South Carolina, USA. Oceanography, 26(3), 78-84.
#' @references Rochlin, I., & Morris, J. T. (2017). Regulation of salt marsh mosquito populations by the 18.6‐yr lunar‐nodal cycle. Ecology, 98(8), 2059-2068.
#' 
predictLunarNodalCycle <- function(year, floodElv, MSL0, MSL, lunarNodalAmp) {
  # Build MHW lines based on MSL, long-term tidal range and lunar nodal amplitude
  MHW <- MSL + (floodElv-MSL0) + (lunarNodalAmp * (sin(2*pi*(year-1983)/18.61)))
  return(MHW)
}