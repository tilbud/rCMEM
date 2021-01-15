## Figure summarising rCMEM behavior. 

require(rCTM)
require(tidyverse)
require(gridExtra)

# Generic input will run from mean sea level to highest tide + sea-level rise 
transectoutput <- runCohortMemTransect(startYear=2015, relSeaLevelRiseInit=0.3, relSeaLevelRiseTotal=100,
                                       initElv=21.9, meanSeaLevel=7.4, 
                                       meanHighWaterDatum=16.9, meanHighHighWaterDatum=25.4,
                                       meanHighHighWaterSpringDatum=31.2, 
                                       suspendedSediment=3e-05, lunarNodalAmp=2.5, bMax=0.25, 
                                       zVegMin=-24.7, zVegMax=44.4, zVegPeak=22.1,
                                       plantElevationType="orthometric", rootToShoot=2,
                                       rootTurnover=0.5, rootDepthMax=30, omDecayRate=0.8,
                                       recalcitrantFrac=0.2, captureRate=2.8,
                                       abovegroundTurnover=1.5,
                                       initElvMin=7,
                                       initElvMax=132,
                                       elvIntervals=4)

scenarioTransectGraph <- transectoutput[[1]] %>% 
  filter(surfaceElevation != initElv) %>% 
  mutate(above_below_msl = ifelse(surfaceElevation >= meanSeaLevel, "above MSL", "below MSL"))

# Which members to visualize in Cohort Graph?

# initElv == 27
# year %in% c(2015, 2065, 2114)
cohortSubset <- transectoutput[[2]] %>% 
  filter(initElv == 27,
         year %in% c(2016, 2065, 2114)) %>% 
  mutate(labs1 = recode(as.character(year),
                        "2016" = "B", 
                        "2065"="C",
                        "2114"="D"),
         labs2 = paste(labs1, as.character(year), sep=". Cohort depth series, "))

scenarioSubset <- transectoutput[[1]] %>% 
  filter(initElv == 27,
         year %in% c(2016, 2065, 2114)) %>% 
  mutate(labs1 = recode(as.character(year),
                        "2016" = "B", 
                        "2065"="C",
                        "2114"="D"),
         labs2 = paste(labs1, as.character(year), sep=". Cohort depth series, "),
         above_below_msl = ifelse(surfaceElevation >= meanSeaLevel, "above MSL", "below MSL")) 
  

scenarioTransect <- ggplot(data=scenarioTransectGraph, aes(x=year, y=surfaceElevation, color=above_below_msl)) +
  geom_line(aes(group=as.character(initElv))) +
  geom_point(aes(shape=above_below_msl)) +
  scale_shape_manual(values=c(24, 25)) +
  geom_point(data = scenarioSubset, size = 1.25, color = "black", aes(shape=above_below_msl)) +
  geom_text(data = scenarioSubset, aes(label = labs1), nudge_x = -2, color = "black") +
  theme_minimal() +
  ylab("Surface Elevation (cm NAVD88)") +
  theme(legend.title = element_blank(), 
        legend.position = "bottom") +
  ggtitle("A. cMEM Behavior Over an Elevation Transect")

(scenarioTransect)


# First reshape the mass cohorts so that they're in long form
mass_cohorts <- cohortSubset %>%
  dplyr::select(-cumCohortVol, -respired_OM) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::group_by(year, initElv) %>% 
  dplyr::mutate(cohortIndex = length(age):1) %>% 
  ungroup() %>% 
  tidyr::gather(key = "mass_pool", value = "mass_fraction", 
                -age, -year, -layer_top, -layer_bottom, -cohortIndex, -initElv, -labs1, -labs2, -elvMin, -elvMax) %>%
  dplyr::group_by(year, age, cohortIndex, initElv) %>%
  dplyr::mutate(mass_pool = str_replace(mass_pool, "_", " "),
                mass_pool = factor(mass_pool, 
                                   levels=c("mineral", 
                                            "slow OM", 
                                            "fast OM", 
                                            "root mass"))) %>%
  dplyr::arrange(year, age, mass_pool) %>%
  dplyr::mutate(max_mass = cumsum(mass_fraction),
                min_mass = ifelse(mass_pool==first(mass_pool),0,lag(max_mass)),
                mass_pool = as.character(mass_pool)) %>%
  # Join mass cohorts with scenario table to convert depths to referenced elevations
  dplyr::ungroup() %>%
  dplyr::left_join(scenarioSubset) %>%
  dplyr::mutate(layer_top=surfaceElevation-layer_top, 
                layer_bottom=surfaceElevation-layer_bottom)

# Reshape the scenario table
tides <- scenarioSubset %>%
  # Track any elevation threholds in the animation speciefied.
  # meanSeaLevel and meanHighWater are the defaults
  dplyr::select(year, meanSeaLevel, meanHighWater, labs2) %>%
  tidyr::gather(value="WaterLevel", key="datum", -year, -labs2) %>%
  # dplyr::rename(year=years) %>%
  dplyr::arrange(year) %>%
  dplyr::filter(complete.cases(.)) %>% 
  mutate(datum = recode(datum, 
                             "meanSeaLevel"= "MSL",
                             "meanHighWater"="MHW"))

# get rid of any NA values.               
mass_cohorts_almostAll <- mass_cohorts %>%  dplyr::filter(complete.cases(.))

chPalette = c("#56B4E9", "#999999", "#E69F00", "#009E73")

# gganimate stuff
graph_mass_cohorts <- ggplot2::ggplot(data = mass_cohorts_almostAll, 
                                        aes(xmin = min_mass, xmax = max_mass, 
                                            ymin = layer_top, ymax = layer_bottom
                                        )) +
  ggplot2::geom_rect(aes(fill = mass_pool), color = rgb(0,0,0, alpha = 0.1)) +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(values=chPalette) +
  ggplot2::geom_hline(data=tides, aes(yintercept=WaterLevel, lty=datum), color="blue") +
  ggplot2::ylab("Depth (cm NAVD88)") +
  ggplot2::xlab("Mass Accumulated Per Cohort (g)") +
  facet_wrap(.~labs2) +
  theme(legend.position = "right",
        legend.title = element_blank())

(graph_mass_cohorts)

grid.arrange(scenarioTransect, graph_mass_cohorts, nrow = 2)

cMemFig <- arrangeGrob(scenarioTransect, graph_mass_cohorts)  

ggsave("temp/cMemBehaviorFig.pdf", width=7.25, height=7.25, dpi=300, cMemFig)
ggsave("temp/cMemBehaviorFig.jpg", width=7.25, height=7.25, dpi=300, cMemFig)
