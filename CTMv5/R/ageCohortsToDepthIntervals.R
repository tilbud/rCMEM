ageCohortsToDepthIntervals <- function(inputDF = startingConditionChohorts,
                                       depthMin.x = "depthMin", 
                                       depthMax.x = "depthMax",
                                       depthMin.y = 0:99, 
                                       depthMax.y = 1:100) {
  
  tempDf <- inputDF %>%
    rename(depthMin = depthMin.x,
           depthMax = depthMax.x)
  
  # First define your cohort depths and age column names
  depthsAndAges <-c("depthMin", "depthMax")
  
  # All other columns are response columns
  targetAttributes <- tempDf %>%
    select(names(tempDf)[! names(tempDf) %in% c(depthsAndAges)])
  
  # Then start stepping one age cohort at a time
  for (targetDepth in 1:length(depthMin.y)) {
    
    targetMin <- depthMin.y[targetDepth]
    targetMax <- depthMax.y[targetDepth]
    
    # Apply the above return_overlap funtion to the temporary data frame
    #   to get a series of weights indicating how much overlap there is between
    #   the chort we're iterating through and the age of the depth intervals in the 
    #   table. These should add up to one. Most should be 0.
    depthWeights <- mapply(return_overlap,
                           x1=targetMin, 
                           x2=targetMax, 
                           y1=tempDf$depthMin, 
                           y2=tempDf$depthMax)
    depthWeights <- depthWeights / sum(depthWeights) # make weights sum to 1
    
    # For each target attribute multiply the depth series by the weights, by transforming
    # it to a matrix, getting the colmun sums and then converting back to a dataframe.
    weightedAttributesMatrix <- colSums(as.matrix(targetAttributes * depthWeights))
    weightedAttributesDf <- as.data.frame(t(weightedAttributesMatrix))
    
    if (targetDepth == 1) {
      outputDF <- weightedAttributesDf
    } else {
      outputDF <- rbind(outputDF, weightedAttributesDf)
    }
  }
  return(outputDF)
}
