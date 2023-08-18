#' Run Multiple Species Biomass Functions
#' 
#' This function takes elevation and biological inputs, runs parabolic predictBiomass function for multiple species, and returns a single above ground biomass value and set of biological parameters based on a competition function.
#' @param z a numeric, elevation
#' @param bMax a numeric, or vector of numerics, maximum biomass
#' @param zVegMax a numeric, or vector of numerics, upper elevation of biomass limit
#' @param zVegMin a numeric, or vector of numerics, lower elevation of biomass limit
#' @param zVegPeak (optional) a numeric, or vector of numerics, elevation of peak biomass
#' @param rootToShoot a numeric, or vector of numerics, root to shoot ratio
#' @param rootTurnover a numeric, or vector of numerics, belowground biomass annual turnover rate
#' @param recalcitrantFrac
#' @param abovegroundTurnover (optional) a numeric, or vector of numerics, aboveground biomass annual turnover rate
#' @param abovegroundSlowpoolFrac
#' @param rootDepthMax a numeric, or vector of numerics, maximum (95\%) rooting depth
#' @param speciesCode (optional) a character, or vector of characters, species names or codes associated with biological inputs
#' @param competition.fn (optional) a function that takes at least a dataframe bio_table as an input, models competition between multiple species, and outputs an one row data frame aggregating biomass and biological inputs
#' 
#' @return a one dataframe with above ground biomass, and biological parameters representing the dominant specie(s) at the elevation
#' @export
runMultiSpeciesBiomass <- function(z, bMax, zVegMax, zVegMin, zVegPeak=NA,
                                   rootToShoot, rootTurnover, recalcitrantFrac, rootDepthMax, 
                                   abovegroundTurnover=1, abovegroundSlowpoolFrac,
                                   speciesCode=NA, competition.fn=NA) {
   # If a custom competition function is not specified ...
   if (is.na(competition.fn)) {
     # Generic competition function filters the inputed bio_table maximum aboveground biomass
     competition.fn <- function(bio_table) {
       return(dplyr::filter(bio_table, aboveground_biomass == max(bio_table$aboveground_biomass)))
     }
   }
  
   # make a list of all the biological inputs
   bio_inputs <- list(bMax, zVegMax, zVegMin, zVegPeak, rootToShoot, 
                     rootTurnover, recalcitrantFrac, rootDepthMax, 
                     abovegroundTurnover, abovegroundSlowpoolFrac, speciesCode)
   
   # Get the length of each input
   input_lengths <- sapply(bio_inputs, length)
   
   # All the lengths either need to be 1 or equal to the number of species that are inputted
   if (all(input_lengths == 1 | input_lengths==max(input_lengths))) {
     # If species codes are not specified ... 
     if (all(is.na(speciesCode)) | length(speciesCode)==max(input_lengths)) {
       if (any(is.na(speciesCode))) {
         # ... then assign them 1,2,3,etc.
         speciesCode <- as.character(1:max(input_lengths))
       }
       # Create a table with all the biological inputs
       bio_table <- data.frame(speciesCode=speciesCode, 
                               stringsAsFactors = F) %>%
         # Using mutate will add vectors for all mutliple values and repeat single values
         # so all the columns will be the same length.
         dplyr::mutate(bMax = bMax,
                       zVegMax=zVegMax,
                       zVegPeak=zVegPeak,
                       zVegMin=zVegMin,
                       rootToShoot=rootToShoot, 
                       rootTurnover=rootTurnover,
                       recalcitrantFrac=recalcitrantFrac,
                       abovegroundTurnover=abovegroundTurnover,
                       abovegroundSlowpoolFrac=abovegroundSlowpoolFrac,
                       rootDepthMax=rootDepthMax)
     } else {
       stop("Species codes either need to be blank, or have the same number as biomass inputs.")
     }
   } else {
     stop("The number of biomass inputs either need to be the same length for multiple species or singular.")
   }
   
   # Run the parabolic biomass function on the table
   bio_table$aboveground_biomass <- mapply(predictBiomass, z=z, 
                                          bMax=bio_table$bMax,
                                          zVegMax = bio_table$zVegMax, 
                                          zVegMin = bio_table$zVegMin, 
                                          zVegPeak = bio_table$zVegPeak)
   
   bio_table$belowground_biomass <- bio_table$aboveground_biomass * bio_table$rootToShoot
   
   # If all aboveground biomass values are 0 ...
   if (all(bio_table$aboveground_biomass == 0)) {
     # ... then make all bio params 0 and changes species name to unvegetated.
     bio_table[,6:11] <- 0
     bio_table[,1] <- "unvegetated"
     bio_table <- bio_table[1,]
   }
   
   # If the returned dataframe is more than one row long ...
   if (nrow(bio_table)>1) {
     
     # ... then run the competition function.
     bio_table <- competition.fn(bio_table)
  
     }
     
   # If more than one agb have exactly the same value ... 
   if (nrow(bio_table)>1) {
     # ... then simplify the table to being one row.
     bio_table <- bio_table %>% 
       dplyr::group_by() %>% 
       # Group the species name into a single string ...
       dplyr::summarise(speciesCode = paste(bio_table$speciesCode, collapse="; "), 
                        # ... and average all the parameters.
                        rootToShoot=mean(rootToShoot), 
                        rootTurnover=mean(rootTurnover),
                        recalcitrantFrac=mean(recalcitrantFrac),
                        abovegroundTurnover=mean(abovegroundTurnover),
                        abovegroundSlowpoolFrac=mean(abovegroundSlowpoolFrac),
                        rootDepthMax=mean(rootDepthMax),
                        aboveground_biomass=first(aboveground_biomass),
                        belowground_biomass=mean(belowground_biomass)) %>% 
       dplyr::ungroup()
   }
   
   bio_table <- bio_table[, ! names(bio_table) %in% c("bMax","zVegMax","zVegPeak","zVegMin")]

   # Return the data frame.
   return(bio_table)
}
