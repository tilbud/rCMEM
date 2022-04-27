# R Cohort Marsh Equilibrium Model

Tidal marshes maintain some capacity to gain elevation and to sequester carbon in the face of accelerating sea-level rise. These concepts have been conceptualized over time into the Marsh Equilibrium Model (MEM), which describes the biophysical relationships between plant production, inundation, and marsh elevation, and the Cohort Theory Model (CTM) which tracks organic and inorganic mass pools in age-depth cohorts below the marsh surface. Here we present an open source numerical versions of a combined model: the R Cohort Marsh Equilibrium Model \pkg{rCMEM}.The package contains tools for hindcasting and forecasting tidal marsh elevation, soil structure, and soil carbon flux changes in response to sea-level rise as well as tools for visualizing model outputs as animations and translating model them to compare predictions to real life data.

## Example to Run rCMEM

```
# Run rCMEM example
cohortMemExample <- runCohortMem(startYear=2000, 
                                 relSeaLevelRiseInit=0.33, # cm per year
                                 relSeaLevelRiseTotal=71, # cm total
                                 initElv=50.2, # cm NAVD88
                                 meanSeaLevel=-6.9, # cm NAVD88
                                 meanSeaLevelDatum =-10.3, # cm NAVD88
                                 meanHighWaterDatum=58.4, # cm NAVD88
                                 meanHighHighWaterDatum=79.0, # cm NAVD88
                                 meanHighHighWaterSpringDatum=99.2, # cm NAVD88
                                 suspendedSediment=3e-05, # grams per cubic centimeter
                                 lunarNodalAmp=2.24, # cm
                                 lunarNodalPhase=1.46, # year
                                 bMax=0.0867, # grams per square centimeter
                                 zVegMin=-0.47, # dimensionless
                                 zVegMax=2.08, # dimensionless
                                 zVegPeak=0.83, # dimensionless
                                 plantElevationType="dimensionless", # specifies dimensionless
                                 rootToShoot=2, # grams per gram
                                 rootTurnover=0.5, # per year
                                 rootDepthMax=30, # cm
                                 omDecayRate=0.5, # per year
                                 recalcitrantFrac=0.2, # fraction
                                 captureRate=2.8 # per tidal cycle
                                 )

# look at the structure of the function output
str(cohortMemExample)

```

## To Uninstall and Reinstall the Latest Version
```
# Uninstall and reinstall latest branch from GitHub

# 1. If rCMEM is loaded and in the memory, forget rCMEM
if ("rCMEM" %in% (.packages())){
  detach("package:rCMEM", unload=TRUE) 
}

# 2. If remotes is not already installed, install it
if (! ("remotes" %in% installed.packages())) {
  install.packages("remotes")
}

# 3. Install package from developer branch of GitHub
devtools::install_github("https://github.com/tilbud/rCMEM")

# 4. Load version into memory
library(rCMEM)

```
## Documentation

For documentation please refer to an extensive [supplemental vignette](/vignettes/CMEM_Vignette_Supplement.Rmd).

