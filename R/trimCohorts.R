#' Trim cohorts
#'
#' Function cleans up no data values from cohorts table
#' @param cohorts data frame tracking soil cohort mass pools
#' @export
trimCohorts <- function(cohorts) {
  cohorts <- cohorts %>%
    dplyr::filter(cumCohortVol!=0)
  return(cohorts)
}