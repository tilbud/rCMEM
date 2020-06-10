#' Deliver sediment flexibly
#'
#' Function decides which sediment module to use based on available inputs
#' @param z a numeric, elevation of the marsh
#' @param ssc a numeric, suspended sediment concentration of the water column
#' @param MSL a numeric, Mean Sea Level
#' @param MHW a numeric, Mean High Water level
#' @param MHHW a numeric, Mean Higher High Water level
#' @param MHHWS a numeric, Mean Higher High Spring Tide Water level
#' @param settlingVelocity a numeric, the number of times a water column will clear per tidal cycle
#' @param ...
#' 
#' @return a numeric, sediment delivered in a year of flooding
#' @export
deliverSedimentFlexibly <- function(z, ssc, MSL, MHW, MHHW=NA, MHHWS=NA, settlingVelocity, ...) {
  # If the user does not include detailed tidal datums
  deliveredSSC <- ifelse(any(is.na(c(MHHW, MHHWS))),
                         # run the simple SSC module
                         deliveredSedimentSimple(z=z, 
                                                 ssc=ssc, 
                                                 MSL=MSL, MHW=MHW, 
                                                 settlingVelocity=settlingVelocity),
                         # If they do include them, run the 3 tide stage module
                         deliveredSediment3TidalCycle(z=z, 
                                                      ssc=ssc, 
                                                      MSL=MSL, 
                                                      MHW=MHW, 
                                                      settlingVelocity=settlingVelocity,
                                                      MHHW=MHHW, 
                                                      MHHWS=MHHWS))
  
  return(deliveredSSC)
  
}