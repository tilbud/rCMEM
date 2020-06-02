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

