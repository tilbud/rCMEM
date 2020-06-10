depthIntervalsToAgeCohorts <- function(inputDF,
                                       ageColumn = "age",
                                       minColumn = "layer_top",
                                       maxColumn = "layer_bottom") {
  # This function takes a table with an age depth model as an attribute,
  # depthMin and depthMax need to be already present.
  # The function iterates through the age depth model to create an annual cohort
  # table with depthMin, depthMax and weighted averages of the rest of the attributes.
  
  # First create a new data frame with min and max ages
  tempDf <- inputDF %>%
    rename(ageMax = ageColumn, # a lot of this could be cleaned up with better data management
           depthMin = minColumn,
           depthMax = maxColumn) %>%
    mutate(ageMin = ifelse(ageMax == min(ageMax),
                           0, lag(ageMax)),
           accretionRate = (depthMax - depthMin)/(ageMax-ageMin)) %>%
    select(depthMin, depthMax, ageMin, ageMax, everything())
  
  # Define cohort depths and age column names
  depthsAndAges <-c("depthMin", "depthMax", "ageMin", "ageMax")
  
  # All other columns are targeted for weighted averaging.
  targetAttributes <- tempDf %>%
    select(names(tempDf)[! names(tempDf) %in% c(depthsAndAges)])
  
  # Next, determine how many cohorts you need. Round up.
  nCohorts <- ceiling(max(inputDF[ageColumn]))
  
  # Then start stepping one age cohort at a time
  for (cohortN in 1:nCohorts) {
    
    # Define the cohort's start and end year
    cohortStartYear <- cohortN-1
    cohortEndYear <- cohortN
    
    # Apply the above return_overlap funtion to the temporary data frame
    #   to get a series of weights indicating how much overlap there is between
    #   the cohort we're iterating through and the age of the depth intervals in the 
    #   table. These should add up to one. Most should be 0.
    depthWeights <- mapply(return_overlap,
                           x1=cohortStartYear, x2=cohortEndYear, 
                           y1=tempDf$ageMin, y2=tempDf$ageMax)/(tempDf$ageMax-tempDf$depthMin)
    
    
    # For each target attribute multiply the depth series by the weights, by transforming
    # it to a matrix, getting the colmun sums and then converting back to a dataframe.
    weightedAttributesMatrix <- colSums(as.matrix(targetAttributes * depthWeights))
    weightedAttributesDf <- as.data.frame(t(weightedAttributesMatrix))
    
    # If first cohort, create an output data frame; else add temporary data to 
    # output using a row bind.
    if (cohortN == 1) {
      outputDF <- weightedAttributesDf
      cohort_depth_max <-c(weightedAttributesDf$accretionRate) 
    } else {
      outputDF <- rbind(outputDF, weightedAttributesDf)
      cohort_depth_max <-c(cohort_depth_max,
                           cohort_depth_max[cohortN-1] + weightedAttributesDf$accretionRate) 
    }
  }
  # Clean up output table and make sure all needed columns are present and in order.
  outputDF <- outputDF %>%
    mutate(cohort = nCohorts:1,
           age_from_surface = 1:nCohorts,
           depthMax = cohort_depth_max,
           depthMin = ifelse(age_from_surface == 1, 0, lag(depthMax))
    ) %>%
    select(cohort, age_from_surface, depthMin, depthMax, everything())
  return(outputDF)
}
