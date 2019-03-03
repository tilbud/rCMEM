#' Fill depth increments using simple Cohort Theory Model.
#' 
#' This function fills 1 cm depth increments down to a maximum depth interval 
#' using a simple version of Cohort Theory Model. It takes several biotic and 
#' abiotic inputs, and fills cells assuming steady state conditions.
#' It outputs a table of depth increments with root mass, refractory 
#' and labile organic matter, inorganic matter, and age-depth relationships.
#' There are optionalities to create derivative statistics used to simulate
#' measurements and to graph the outputs.

#' @param totalDepth 
#' @param rootDepthMax Depth of the root zone in cm below surface.
#' @param totalRootBmass Total Root Biomass, g/m2.
#' @param refractoryFrac Refractory Portion of organic production, g/g.
#' @param bgTurnover  Below Ground Turnover Rate, 1/yr.
#' @param ssc Suspended Sediment Concentration, mg per liter.
#' @param depthBelowMHW Depth of Marsh Surface Below Mean High Water.
#' @param omDecayRate Labile organic matter decay rate, 1/yr.
#' @param coreYear Core Collection Year CE.
#' @param rootShape Could be linear, exponential, or custom.
#' @param k1 Organic and Inorganic Self Packing Densities: k1 and k2.
#' @param k2 
#' @param extraStatsOn 
#' @param graphingOn 
#' @param plotName 
#'
#' @return
#' @export
#'
#' @examples
fillDepthsUsingSimpleCTM <- function(
  totalDepth = 150,
  rootDepthMax = 30, # Depth of the root zone in cm below surface
  totalRootBmass = 3000, # Total Root Biomass, g/m2
  refractoryFrac = 0.2, # Refractory Portion of organic production, g/g
  bgTurnover = 0.5,  # Below Ground Turnover Rate, 1/yr
  ssc = 20, # Suspended Sediment Concentration, mg per liter
  depthBelowMHW = 10, # Depth of Marsh Surface Below Mean High Water
  omDecayRate = 0.8, # Labile organic matter decay rate, 1/yr
  coreYear = 2015, # Core Collection Year CE
  rootShape = "linear", # Could be linear, exponential, or custom
  k1 = 0.085,   # Organic and Inorganic Self Packing Densities: k1 and k2
  k2 = 1.99,
  extraStatsOn = F,
  graphingOn = F,
  plotName = "rCtmOutput") {
  
  # Conversions and constants
  nTidesPerYear <- 704
  meanTidalHeight <- depthBelowMHW
  ssc_gPerCm2 <- ssc * 0.000001 # convert mg/l to grams/cm^2
  cumAnnWaterVol <- nTidesPerYear * meanTidalHeight # Cumulative water volume
  annSediment <- ssc_gPerCm2 * cumAnnWaterVol
  Csyr <- coreYear - 1963 # Year relative to peak Cesium-137 deposition
  
  # Internal function for calculating years of input at steady state 
  # given a root mass and an annual suspended sediment deposition.
  calculateYearsOfInput <- function(x) {
    # For the solver the function needs to have an input x and output y.
    yrsOfInput <- x
    
    mineralSection <- yrsOfInput * annSediment
    
    # If this is the first depth cohort
    if (section == 1) {
      # Annual refractory productivity in section.
      refracOM_temp <- yrsOfInput * refractoryFrac * bgTurnover * rootsInSection
      
      # Annual labile productivity in section
      Labprod_temp = (1 - refractoryFrac) * bgTurnover * rootsInSection
      
      # Original labile OM remaining at time t
      LabileOM_temp = yrsOfInput * Labprod_temp * (1 - omDecayRate) ^ yrsOfInput
      
    } else {
      # If this is an older/deeper cohort at steady state carry over
      # refractory organic matter production from the previous  cohort.
      refracOM_temp = yrsOfInput * refractoryFrac * bgTurnover * 
        rootsInSection + refracOM[section - 1]
      
      # Annual labile production,
      Labprod_temp = (1 - refractoryFrac) * bgTurnover * rootsInSection
      
      # If this is an older/deeper cohort at steady state carry over
      # labile organic matter productuvity from previous cohort.
      LabileOM_temp = yrsOfInput * Labprod_temp * 
        (1 - omDecayRate) ^ yrsOfInput + 
        LabileOM[section - 1] * (1 - omDecayRate) ^ yrsOfInput
    }
    
    # Calculate section volume given years of accumulation
    # using ideal mixing model.
    cohortVol <- mineralSection / k2 + 
      (rootsInSection + LabileOM_temp + refracOM_temp) / k1
    
    y = cohortVol - 1 # set volume equal to 0
    
    # Return 
    return(y)
  }
  
  # Create output table and fills cells with roots
  outputTable <- fillDepthCellsWithRoots(depthMins = 0:(totalDepth-1),
                                         depthMaxs = 1:totalDepth,
                                         totalRootBmass = totalRootBmass,
                                         rootShape = rootShape,
                                         rootDepthMax = rootDepthMax
  )
  # Convert to grams per cm^2
  rootBmass_gPerCm2 <- outputTable$rootBiomass_gPerM2 / 10000
  
  # Define a steady state accretion rate using ideal mixing model
  accRate <- (totalRootBmass * omDecayRate * bgTurnover) / 
    k1 + annSediment / k2
  
  # First guess accretion rate using ideal mixing model
  y2start <- 1 / (5 * ((rootBmass_gPerCm2[1] * omDecayRate * bgTurnover) / k1 
                       + annSediment / k2))
  if (y2start == 0) { y2start <- 0.8 / accRate}
  
  # Empty vectors for temporary storage
  refracOM <- c()
  Labprod <- c()
  LabileOM <- c()
  
  # Begin the iterations over n sections
  for (section in 1:nrow(outputTable)) {
    
    # For each iteration add roots per section
    rootsInSection <- rootBmass_gPerCm2[section]
    
    # Use a nonlinear solver to calculate years per depth interval using 1st 
    # guess and predifined function. Ouput is years of input
    yrsOfInput <- nleqslv(y2start, calculateYearsOfInput)$x
    
    # Calulate mineral mass
    mineralSection <- yrsOfInput * annSediment
    
    # If this is the first depth cohort
    if (section == 1) {
      # Annual refractory productivity in section.
      refracOM_temp <- yrsOfInput * refractoryFrac * bgTurnover * rootsInSection
      
      # Annual labile productivity in section
      Labprod_temp = (1 - refractoryFrac) * bgTurnover * rootsInSection
      
      # Original labile OM remaining at time t
      LabileOM_temp = yrsOfInput * Labprod_temp * (1 - omDecayRate) ^ yrsOfInput
      
    } else {
      # If this is an older/deeper cohort at steady state carry over
      # refractory organic matter production from the previous  cohort.
      refracOM_temp = yrsOfInput * refractoryFrac * bgTurnover * 
        rootsInSection + refracOM[section - 1]
      
      # Annual labile production,
      Labprod_temp = (1 - refractoryFrac) * bgTurnover * rootsInSection
      
      # If this is an older/deeper cohort at steady state carry over
      # labile organic matter productuvity from previous cohort.
      LabileOM_temp = yrsOfInput * Labprod_temp * 
        (1 - omDecayRate) ^ yrsOfInput + 
        LabileOM[section - 1] * (1 - omDecayRate) ^ yrsOfInput
    }
    
    # Add to memory
    refracOM[section] <- refracOM_temp
    Labprod[section] <- Labprod_temp
    LabileOM[section] <- LabileOM_temp
    
    # Update output table
    outputTable$yearsOfInput[section] <- yrsOfInput
    outputTable$mineral_gPerCm2[section] <- mineralSection
    
    # Iterate through next section
    
  }
  # Done adding roots to all depth intervals
  # Update output table
  outputTable$refractoryBiomass_gPerM2 <- 10000 * refracOM
  outputTable$labileOM_gPerM2 <- 10000 * LabileOM
  
  # Use some dplyr operations to get age depth and put columns in order
  outputTable <- outputTable  %>%
    mutate(horizonAge_yrs = cumsum(yearsOfInput)) %>%
    select(depthMin, depthMax, horizonAge_yrs, yearsOfInput, rootBiomass_gPerM2, 
           mineral_gPerCm2, refractoryBiomass_gPerM2, labileOM_gPerM2)
  
  # Create Output Table with Summary Variables
  outputTable2 <- data.frame(variable = rep(NA, 6),
                             value = rep(NA, 6),
                             units = rep(NA, 6))
  
  # Create Table Comparing the Numerical and Analytic Quantities
  outputTable3 <- data.frame(method = c("Numerical", "Analytical"),
                             decay = c(NA, NA),
                             BGProd = c(NA, NA),
                             refractoryOMProd = c(NA, NA))
  
  # If we turned on graphing or extra stats run this extra stats module.
  # This can be turned off to speed up the code.
  if (extraStatsOn == T | graphingOn == T) {
    
    # Use dplyr operations to get 210Pb profiles
    outputTable <- outputTable %>%
      mutate(Pb210_dmpPerG = ifelse(mineral_gPerCm2 > 0,
                                    5 * yearsOfInput * exp(-horizonAge_yrs * 0.69315 / 22.3) / mineral_gPerCm2,
                                    NA),
             lnPb210 = log(Pb210_dmpPerG)
      )
    
    # Use vector subsetting to get 137Cs based marker position and accretion
    secCsyr <- which(outputTable$horizonAge_yrs >= Csyr)[1]
    acrCsyr <- secCsyr / outputTable$horizonAge_yrs[secCsyr]
    
    # Use dplyr operations to calculate loss on ignition and bulk density
    outputTable <- outputTable %>%
      mutate(lossOnIgnition = 0.9 * (rootBiomass_gPerM2 + refractoryBiomass_gPerM2 + labileOM_gPerM2) / 
               (mineral_gPerCm2*10000 + rootBiomass_gPerM2 + refractoryBiomass_gPerM2 + labileOM_gPerM2),
             dryBulkDensity_gPerCc = 1 / (lossOnIgnition / k1 + (1 - lossOnIgnition) / k2))
    
    # Add variables to output table 2 for export and display
    outputTable2$variable[1] <- "Steady State Accretion"
    outputTable2$value[1] <- 10 * (refractoryFrac * bgTurnover * totalRootBmass * 0.0001 / k1 + annSediment / k2)
    outputTable2$units[1] <- "mm/yr"
    
    outputTable2$variable[2] <- "Surface Accretion"
    outputTable2$value[2] <- 1 / outputTable$yearsOfInput[1] * 10
    outputTable2$units[2] <- "mm/yr"
    
    # Calculate Pb210 Accretion
    PbModel <- lm(x~y, data = data.frame(
      x = outputTable$lnPb210[1:80], y = outputTable$depthMax[1:80]))
    PbSlope <- PbModel[[1]]["y"]
    
    Dhalf = -0.69315 / PbSlope
    PBAR = Dhalf / 22.3
    
    outputTable2$variable[3] <- "Accretion at Cs pk"
    outputTable2$value[3] <- acrCsyr * 10
    outputTable2$units[3] <- "mm/yr"
    
    outputTable2$variable[4] <- "Pb210 Accretion Rate"
    outputTable2$value[4] <- PBAR * 10
    outputTable2$units[4] <- "mm/yr"
    
    # Calculate the weighted organic matter age at the cs peak
    yrsold <- c()
    for (i in 1:secCsyr) {
      yrsold[i] <- 0
      for (j in i:secCsyr) {
        yrsold[i] = yrsold[i] + outputTable$yearsOfInput[j]
      }
      
      organicMatterAge <- 0
      for (j in 1:secCsyr) {
        organicMatterAge <- organicMatterAge + 
          (yrsold[j] * outputTable$yearsOfInput[j] * rootBmass_gPerCm2[j] * refractoryFrac * bgTurnover)
      }
      
      wtage <- (organicMatterAge + 
                  (outputTable$yearsOfInput[secCsyr] * rootBmass_gPerCm2[secCsyr]) + 
                  (outputTable$yearsOfInput[secCsyr] * LabileOM[secCsyr])) / 
        (refracOM[secCsyr] + LabileOM[secCsyr] + rootBmass_gPerCm2[secCsyr])
    }
    
    # Add outputs to output tables.
    outputTable2$variable[5] <- "Weighted C Age"
    outputTable2$value[5] <- wtage
    outputTable2$units[5] <- "yr"
    
    outputTable2$variable[6] <- "Cesium Horizon Age"
    outputTable2$value[6] <-  round(outputTable$horizonAge_yrs[secCsyr])
    outputTable2$units[6] <- paste("age @ ", secCsyr, " cm", sep="")
    
    # Make outputs round numbers for display purposes.
    outputTable2 <- outputTable2 %>%
      mutate_if(is.numeric, round, 1)
    
    # Compare numerical vs analytics solutions.
    # First use dplyr operations to calculate decay per section.
    outputTable <- outputTable %>%
      mutate(decayPerSection_gPerM2PerYr = ifelse(depthMax == 1,
                                                  Labprod * omDecayRate,
                                                  Labprod * omDecayRate + lag(LabileOM) * omDecayRate),
             decayPerSection_gPerM2PerYr = decayPerSection_gPerM2PerYr * 10000
      )
    
    # Numerical annual decay 
    outputTable3$decay[1] <- sum(outputTable$decayPerSection_gPerM2PerYr) 
    
    # Analytical Annual decay
    outputTable3$decay[2] <- (1 - refractoryFrac) * bgTurnover * totalRootBmass
    
    # Numerical Below Ground Productivity
    outputTable3$BGProd[1] <- sum(outputTable$rootBiomass_gPerM2 * bgTurnover)
    
    # Analytical Below Ground Productivity
    outputTable3$BGProd[2] <- totalRootBmass * bgTurnover
    
    # Note the following line was corrected from Jim's script
    # Numerical Refractory Organic matter Productivity
    outputTable3$refractoryOMProd[1] <- sum(outputTable$rootBiomass_gPerM2 * bgTurnover * refractoryFrac)
    
    # Analytical Organic matter Productivity
    outputTable3$refractoryOMProd[2] <- totalRootBmass * bgTurnover * refractoryFrac
    
    # Round digits for display purposes
    outputTable3 <- outputTable3 %>%
      mutate_if(is.numeric, round)
    
    # If graphing is on, create plots. This can be turned off to speed up code.
    if (graphingOn == T)  {
      
      # Create nice display text tables.
      inputTable <- as_tibble(data.frame(
        inputName = c("rootDepthMax", "totalRootBmass", "refractoryFrac", "ssc", 
                      "omDecayRate", "bgTurnover", "depthBelowMHW", "coreYear", 
                      "k1", "k2"),
        inputValue = c(rootDepthMax, totalRootBmass, refractoryFrac, ssc, 
                       omDecayRate, bgTurnover, depthBelowMHW, coreYear, k1, k2),
        units = c("cm","g/m2","g/g","mg/liter", "1/yr","1/yr","cm3","date",
                  "g/cm3","g/cm3")
      ))
      
      displayTable1 <- ggtexttable(inputTable, rows = NULL) 
      
      displayTable2 <- ggtexttable(outputTable2, rows = NULL) 
      
      displayTable3 <- ggtexttable(outputTable3, rows = NULL) 
      
      # Create figures using ggplot.
      graph4 <- ggplot(data = outputTable, aes(x = depthMax, y = Pb210_dmpPerG)) +
        geom_line(col = "darkgrey") +
        geom_point(pch = 21, bg = "darkgrey", color = "red") +
        xlim(0, 80) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab(expression(paste({}^"210", "Pb (dpm/g)", sep=""))) +
        scale_y_log10(limits=c(0.1, 100)) +
        ggtitle(expression(paste("Excess "^"210", " Pb (dpm/g)", sep="")))
      
      graph5 <- ggplot(data = outputTable, aes(x = depthMax, y = horizonAge_yrs)) +
        geom_line(col = "darkgrey") +
        geom_point(pch = 21, bg = "darkgrey", color = "red") +
        xlim(0, 80) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab("Horizon Age (years)") +
        ggtitle("Horizon Age")
      
      graph6 <- ggplot(data = outputTable, aes(x = depthMax, y = lossOnIgnition)) +
        geom_line(col = "darkgrey") +
        geom_point(pch = 21, bg = "darkgrey", color = "red") +
        xlim(0, 50) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab("OM Concentration (%)") +
        ggtitle("OM Concentration (% of Dry Weight)")
      
      graph7 <- ggplot(data = outputTable, aes(x = depthMax, y = mineral_gPerCm2)) +
        geom_line(col = "darkgrey") +
        geom_point(pch = 21, bg = "darkgrey", color = "red") +
        xlim(0, 50) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab(expression(paste("Inorganic Sediment (g/cm"^"2",")", sep=""))) +
        ggtitle(expression(paste("Mineral Content", sep="")))
      
      graph8 <- ggplot(data = outputTable, aes(x = depthMax, y = rootBiomass_gPerM2)) +
        geom_line(col = "darkgrey") +
        geom_point(pch = 21, bg = "darkgrey", color = "red") +
        xlim(0, 50) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab(expression(paste("Biomass (g/cm"^"2",")", sep=""))) +
        ggtitle(expression(paste("Live Root Biomass", sep="")))
      
      carbon_table <- outputTable %>%
        select(depthMax, labileOM_gPerM2, refractoryBiomass_gPerM2) %>%
        rename(labile = labileOM_gPerM2,
               refractory = refractoryBiomass_gPerM2) %>%
        gather(key = OM_Status, value = gPerM2, labile, refractory)
      
      graph9 <- ggplot(data = carbon_table, aes(x = depthMax, y = gPerM2)) +
        geom_line(aes(color = OM_Status)) +
        geom_point(pch = 21, bg = "darkgrey", aes(color = OM_Status)) +
        xlim(0, 50) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab(expression(paste("Organic Matter (g/cm"^"3",")", sep=""))) +
        labs(color = "") +
        theme(legend.position = "top")
      
      graph10 <- ggplot(data = outputTable, aes(x = depthMax, y = dryBulkDensity_gPerCc)) +
        geom_line(col = "darkgrey") +
        geom_point(pch = 21, bg = "darkgrey", color = "red") +
        xlim(0, 50) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab(expression(paste("Bulk Density (g/cm"^"3",")", sep=""))) +
        ggtitle(expression(paste("Dry Bulk Density", sep="")))
      
      graph11 <- ggplot(data = outputTable, aes(x = depthMax, y = decayPerSection_gPerM2PerYr)) +
        geom_line(col = "darkgrey") +
        geom_point(pch = 21, bg = "darkgrey", color = "red") +
        xlim(0, 50) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab(expression(paste("Decay (g m"^"-2", " yr"^"-1",")", sep=""))) +
        ggtitle(expression(paste("Annual Decay", sep="")))
      
      graph12 <- ggplot(data = outputTable, aes(x = depthMax, y = yearsOfInput)) +
        geom_line(col = "darkgrey") +
        geom_point(pch = 21, bg = "darkgrey", color = "red") +
        xlim(0, 90) +
        theme_minimal() +
        xlab("Depth (cm)") +
        ylab("Years") +
        ggtitle("Years of Input per 1 cm section\nor # of cohorts per section")
      
      # Set output pdf location, name, and dimensions.
      jpeg(paste("temp/", plotName, ".jpg", sep=""), width=11, height=8.5, units="in", res=300)
      
      # Arrange all tables and figures into a single plot and print it to pdf.
      outputFullTable <- grid.arrange(nrow = 4, 
                                      displayTable1, displayTable3, displayTable2,
                                      graph4, graph5, graph6,
                                      graph7, graph8, graph9,
                                      graph10, graph11, graph12,
                                      heights = c(1.75,1,1,1))
      (outputFullTable)
      dev.off() # Turn graphics off
    }
  }
  return(list(outputTable, outputTable2, outputTable3))
}
