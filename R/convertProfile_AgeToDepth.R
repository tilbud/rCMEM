#' Convert soil profile from age to depth based cohorts
#' 
#' This function takes an age cohort soil profile and converts it to a depth based cohort using weighted sums and means.
#'
#' @param ageCohort a dataframe specifying the age cohort. Must contain \code{layer_top}, \code{layer_bottom}, \code{age}, \code{fast_OM}, \code{slow_OM}, \code{mineral}, and \code{root_mass}.
#' @param layerTop a vector of the tops of the soil layers we want to convert to.
#' @param layerBottom a vector of the bottom of the soil layers we want to convert to.
#'
#' @return a data frame specifying the depth cohoprts. Contains \code{layer_top}, \code{layer_bottom}, \code{age}, \code{fast_OM}, \code{slow_OM}, \code{mineral}, and \code{root_mass}.
#' @export
#'
#' @importFrom plyr ddply
convertProfile_AgeToDepth <- function(ageCohort, layerTop, layerBottom){
  
  ans <- plyr::ddply(data.frame(bottom=layerBottom, top=layerTop), c('top', 'bottom'), function(xx){
    
    layerWeights <- pmax( pmin(ageCohort$layer_bottom, xx$bottom) - pmax(ageCohort$layer_top, xx$top), 
                          0) / (ageCohort$layer_bottom - ageCohort$layer_top)
    
    ans <- base::lapply(ageCohort[,c('fast_OM', 'slow_OM', 'mineral', 'root_mass')],
                        function(yy)sum(yy*layerWeights))
    
    ans$age <- weighted.mean(ageCohort$age, layerWeights)
    ans$input_yrs <- sum(layerWeights)
    ans$layer_bottom <- xx$bottom
    ans$layer_top <- xx$top
    
    return(as.data.frame(ans))
  })
  
  return(ans[,c('layer_top', 'layer_bottom', 'age', 'input_yrs',
                'fast_OM', 'slow_OM', 'mineral', 'root_mass')])
}