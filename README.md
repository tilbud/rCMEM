# Cohort Theory and Marsh Equilibrium Models

This repository corresponds to the R package in development for use of the Cohort Theory Model (CTM) which is a simple method for accounting for volume and mass changes in tidal systems. This package also contains a special implementation of CTM which uses the Marsh Equilibrium Model (MEM) to drive inorganic surface sedimentation and below ground organic matter addition, dynamically as a function of sea-level.

# To Download This Developer's Branch

## To Uninstall and Reinstall the Latest Version
```
# Uninstall and reinstall developer branch from GitHub

# 1. If rCTM is already installed, uninstall rCTM
if ("rCTM" %in% installed.packages()) {
  remove.packages('rCTM')
}

# 2. If rCTM is loaded and in the memory, forget rCTM
if ("rCTM" %in% (.packages())){
  detach("package:rCTM", unload=TRUE) 
}

# 3a. If devtools is not already installed, install it
if (! ("devtools" %in% installed.packages())) {
  install.packages('devtools')
}

# 3. Install package from developer branch of GitHub
devtools::install_github("https://github.com/tilbud/rCTM/tree/JimH-dev")

# 4. Load most current version in the memory
library(rCTM)

```

## To Make Sure You Have The Needed Dependencies
```
# Check and install dependencies

# A column vector listing the packages needed
dependencies <- c("tidyverse", 
                  "gganimate", 
                  "gifski",
                  "png")

# Iterate through the necessary packages
for (i in 1:length(dependencies)) {
  print(paste(dependencies[i], " ...", sep=""))
  # If the package is not in the installed package ...
  if (! (dependencies[i] %in% installed.packages())) {
    # .. install it.
    install.packages(dependencies[i])
  } else {
    print("... already installed.")
  }
}

```

## Quick Example to Run MEM with Cohorts and Animate

```
# Run rMEM with Cohorts example
memCohortExample <- runMemWithCohorts(startYear=2015, rslrT1=0.3, rslrTotal=100,
                                      initElv=21.9, MSL=7.4, MHW=16.9, MHHW=25.4, MHHWS=31.2, 
                                      ssc=3e-05, lunarNodalAmp=2.5, bMax=0.25, 
                                      zVegMin=-24.7, zVegMax=44.4, zVegPeak=22.1,
                                      plantElevationType="orthometric", rootToShoot=2,
                                      rootTurnover=0.5, rootDepthMax=30, omDecayRate=0.8,
                                      recalcitrantFrac=0.2, settlingVelocity=2.8,
                                      coreYear = 2050)

# look at the structure of the function output
str(memCohortExample)

# Look at the three tables making up the output
head(memCohortExample$annualTimeSteps)
head(memCohortExample$cohorts)
head(memCohortExample$core)

# run the animate cohorts function
# This will take a few minutes to run and then save an animation of accumulating cohorts
animateCohorts(cohorts=memCohortExample$cohorts,
               scenario=memCohortExample$annualTimeSteps,
               filename="MEM-Try_191212.gif")
```

#File descriptions
1) addCohort
  + Description : Adds another cohort to the soil profile
  + Dependencies : depthOfNonRootVolume, massLiveRoots, sedimentInputs
  + Inputs : massPools, packing densities, time step in years, root turnover time in faction/year, allocation fraction between fast and slow organic matter pool in fraction/year, massLiveRoots, depthOfNonRootVolume, mineral input
  + Outputs : Returns a data frame simular to massPools that includes age, fast_om, slow_OM, mineral, layer top, layer bottoms, and root mass. returns data fram with the added cohort and simulated change in the cohort pools.

2) availableSediment
  + Description : calculates the available sediment for a marsh elevation over a tidal cycle
  + Dependencies : NA
  + Inputs : time per tidal cycle the marsh is inundated, suspended sedument concentration of the water column, Settling velocity - the number of times a water column will clear per tidal cycle, the fraction of available sediment captured by the marsh
  + Outputs : the available sediment given a set of sediment and flooding conditions

3) builHighTideScenario
  + Description : incorporates high tides into a scenario table
  + Dependencies : predictLunarNodalCycle
  + Inputs : mean sea level over tidal datum period, mean high water level, mean higher high water level?,  mean higher high spring tide water level, amplitude of 18 year lunar nodal cycle
  + Outputs : a data frame, including the sea-level and suspended sediment concentraiton scenario inputted, with annual high tide datum(s) added

4) buildScenarioCurve
  + Description : builds an annualied set of MEM inputs
  + Dependencies : NA
  + Inputs :start year of scenario, end year of scenario, mean sea level (a numeric or vector), initial rate of relative sea-level rise, total relative sea-level rise, suspended sediment concentration (a numeric average or vector), 
  + Outputs : a data frame including columns for year, sea-level, and suspended sediment concentration, and rows for each year in the scenario

 5) convertProfile_AgeToDepth
  + Description : Converts soil profile from age to depth based cohorts
  + Dependencies : NA
  + Inputs : a data frame of the age cohort, a vector of the tops of the soil layers we want to convert to, a vector of the bottoms of the soil layers we want to convert to
  + Outputs : a data frame of the depth cohorts

 6) deliverSedimentFlexibly
  + Description : Function decides which sediment module to use based on available inputs
  + Dependencies : "deliveredSediment3TidalCycle, deliveredSedimentSimple"
  + Inputs : marsh elevation, suspended sediment concentration,  mean high water level, mean higher high water level,  mean higher high spring tide water level, number of times a water column will clear per tidal cycle, 
  + Outputs :sediment delivered in a year of flooding (numeric)

 7) deliveredSediment3TidalCycle
  + Description : Calculate delivered sediment, three tidal cycle method
  + Dependencies : availableSediment, floodTimeFromDatum, 
  + Inputs :marsh elevation, suspended sediment concentration, mean sea level, mean high water level, mean higher high water level,  mean higher high spring tide water level, mean low water level, mean lower low water level, mean low lower spring tide water level, number of times a water column will clear per tidal cycle, fraction of available sediment captured by marsh
  + Outputs : sediment delivered over the course of a year
 
 8) deliveredSedimentSimple
  + Description : Calculate delivered sediment, simple method
  + Dependencies : availableSediment
  + Inputs : marsh elevation, fraction of time per tidal cycle marsh is inundated, suspended sediment concentration, mean high water level, mean sea level, mean low water level, number of times water column will clear per tidal cycle, fraction of available sediment captured by marsh, flood events per year
  + Outputs : sediment delivered over the course of a year

 9) depthOfNonRootVolume
  + Description : calculate depth of the given soil volume
  + Dependencies : massLiveRoots.fn
  + Inputs : total root mass per area, max root depth, root density, unit length of soil volume, unit width of soil column, root shape, a numerical array of the cumulative non-rooting volume, relative tolerance accepted, verbose
  + Outputs : depth of the specified non-root volume (numeric)

 10) floodTimeFromDatum
  + Description : Calculate flood time (in hours) from datum
  + Dependencies : NA
  + Inputs :marsh elevation, highest water level of the tidal period, lowest water level of the tidal period
  + Outputs : flood time in hours per tidal period

 11) massLiveRoots
  + Description : function to return the mass of living roots between two layers, given mass per area of roots and a maximum rooting depth
  + Dependencies : NA
  + Inputs : an array of depths from soil top to the bottom of the layer, an array of depths from soil top to the top of the layer, total mass per area of roots, unit length of soil area, unit width of soil area, root shape (linear), maximum root depth, ...
  + Outputs : an array of the associated live root mass for each layer specified
 
12) predictLunarNodalCycle
  + Description : builds a high tide level based on the 18.1 year lunar nodal cycle
  + Dependencies : NA
  + Inputs : mean sea-level over the last complete 18 year datum period, mean sea-level for each scenario, high tide level over the last complete 18 year datum period, amplitude of the 18 year lunar nodal cycle
  + Outputs : MHW
 
 13) predictedBiomass
  + Description : predict aboveground biomass given elevation
  + Dependencies : NA
  + Inputs : marsh elevation, maximum biomass, upper elevation of biomass limit, lower elevation of biomass limit, optional elevation of peak biomass
  + Outputs : aboveground biomass (numeric)
 
 14)runMemWithCohorts
  + Description : Run Marsh Equilibrium Model with Cohorts
  + Dependencies : addCohort, buildHighTideScenario, buildScenarioCurve, convertProfile_AgeToDepth, deliverSedimentFlexibly, predictedBiomass, runToEquilibrium, trimCohorts, zToZstar
  + Inputs : start year, end year, initial rate of relative sea level rise, total relative sea level rise over scenario, initial marsh elevation, initial mean sea level, mean sea level over last datum period, mean high water level over last datum period, (optional) mean higher high water level, (optional) mean higher high spring tide water level, suspended sediment concentration, amplitude of the 18 year lunar nodal cycle, the number of times a water column will clear per tidal cycle, maximum biomass, upper elevation of biomass limit, lower elevation of biomass limit, (optional) elevation of peak biomass, elevation reference of the vegetation growing, root to shoot ratio, ground below biomass annual turnover rate, percentage rooting depth, relationship between depth and root mass, annual fractional mass lost, fraction of organic matter resistant to decay, bulk density of pure organic matter, bulk density of pure mineral matter, bulk density of pure root matter, core year, core depth, core sampling depth minimums, core sampling depth maximums, ...
  + Outputs : a list of data frames, including the annualized summaries, mapped cohorts tracked for every years, and if core year is specified, a core.

 15) runToEquilibrium
  + Description : Run to equilbrium with proscribed surface mineral input
  + Dependencies : addCohort
  + Inputs : minimum run time, maximum run time, relative and absolute tolerance
  + Outputs : a data frame by age cohort with the layer top and bottom depths, age, and mass of the fast organic, slow organic, mineral, and root biomass pools OR a list that is a record of age cohorts

 16) sedimentInputs
  + Description : calculates the annual sediment input into the marsh surface via the following formula: number of tides per tear, times the unit area of interest, times the mean tidal volume above the marsh surface in centimeters cubed, times the suspended sediment concentration, times the conversion factor to convert liter to centimeters cubed. This gives us the mineral mass added per year
  + Dependencies : NA
  + Inputs : suspended sediment concetration, mean tidal height, number of tides per year, soil length and soil width
  + Outputs : the mass of mineral added in one year to the top of the marsh 

17) trimCohorts
  + Description : cleans up no data values from cohorts tables
  + Dependencies : NA
  + Inputs : cohorts
  + Outputs : cohorts

18) zStarToZ
  + Description : convert from Z* to elevation
  + Dependencies : NA
  + Inputs : normalized Z* elevation, mean high water tide, mean sea level
  + Outputs : a numeric describing marsh elevation

19) zToZstar
  + Description : convert elevation to Z*
  + Dependencies : NA
  + Inputs : numeric elevation, mean high water tidde, mean sea level
  + Outputs : elevation normalized to the tidal range
  
  
**Citation:**

_Suggested acknowledgments:_

## Relevant resources and literature  

* See Jim Morris' [MEM web interface for](http://129.252.139.114/model/marsh/mem.asp) inspiration

## Equations and Fundamentals of MEM

Guiding princples:  
*  "Equilibrium" in MEM is a target not a description of the current state.  
*  MEM tries to see how simple you can get the model and still be realistic.  
*  MEM ties elevation of marsh surface to ecosystem productivity.  
*  High elevation, stimulates biomass (to a point), and that causes marsh elevation to raise.  
    
### Ideal Mixing Model

A lot of what underpins current versions of the Marsh Equilibrium Model and the Cohort Theory Model is the ideal mixing model, see [Morris et al. 2016](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015EF000334) and (Holmquist et al 2018](https://www.nature.com/articles/s41598-018-26948-7). This describes the density of a soil as the summation of the soils' relative organic and inorganic fractions, as well as the _self-packing density_ of the organic (k1) and inorganic (k2) matter.

$$BD = {1 \over {{OM\over k1}+ {(1−OM) \over k2}}}$$

If you hold area constant, and you assume you know the input rates of organic and inorganic matter, you can use the ideal mixing model to calculate accretion rate ([Morris et al. 2016](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015EF000334), [Morris 2006](https://www.sciencedirect.com/science/article/pii/S0272771406001776),  [Morris 2007](https://link.springer.com/chapter/10.1007/978-1-4020-6008-3_14)).