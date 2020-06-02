library(rCTM)
library(tidyverse)
library(gridExtra)
library(truncnorm)

# Set up a loop 
iterations <- 50

# Set up priors
# AGB
{
  emin_low <- -10.155 
  mu_emin <- -0.081 
  tau_emin<-3.846
  mu_epeak <- 1.646
  tau_epeak <- 3.846
  mu_emax <- 2.636
  tau_emax <- 3.846 
  emax_high <- 6.968
  mu_bmax <- 2500
  tau_bmax <- 0.001
}

# Below Ground
{
  recalcitrant_alpha	<- 4.36014094 
  recalcitrant_beta <- 10.55162529
  root_depth_alpha <- 4.123560128	
  reoot_depth_beta <- 0.149256165
  decomp_alpha <-	5.621813392	
  decomp_beta <- 4.309662728
  root_shoot_alpha <- 2.868323198	
  root_shoot_beta <- 1.294975575
  root_turnover_alpha <-	2.890021263	
  root_turnover_beta <- 2.598049856
}

# SSC
{
  settelingV_mu <- 2.8 # Calibrated from Jim's site
  # Assume wide uncertainty on this since it's a n of 1
  # s.d. is equivalent to mean
  settelingV_var <- 2.8^2 
  # Gamma distributions are always positive.
  # Moment match to estimate gamma dist parameters from mean and variance
  settelingV_alpha <- settelingV_mu^2/settelingV_var
  settelingV_beta <- settelingV_mu/settelingV_var  
}

# Preallocate membory for data frame
total_iterations <- data.frame(bMax=rep(NA, iterations),
           zVegMin=rep(NA, iterations),
           zVegMax=rep(NA, iterations),
           zVegPeak=rep(NA, iterations),
           rootToShoot=rep(NA, iterations),
           rootTurnover=rep(NA, iterations),
           rootDepthMax=rep(NA, iterations),
           omDecayRate=rep(NA, iterations),
           recalcitrantFrac=rep(NA, iterations),
           settlingVelocity=rep(NA, iterations),
           OMAR = rep(NA, iterations),
           OMSR = rep(NA, iterations))

parameter_names <- c("bMax", "zVegMin", "zVegMax", "zVegPeak",
                     "rootToShoot", "rootTurnover", "rootDepthMax",
                     "omDecayRate", "recalcitrantFrac", 
                     "settlingVelocity")
parameter_types <- c("AGB~Elevation", "AGB~Elevation","AGB~Elevation","AGB~Elevation",
                     "BGB", "BGB", "BGB",
                     "Decay", "Decay",
                     "Sedimentation")

# Set up progress bar
pb <- txtProgressBar(min = 0, max = iterations, style = 3)
for (i in 1:iterations) {
  # Random draws of each prior
  total_iterations$zVegMin[i] <- rtruncnorm(1, a=emin_low, mean=mu_emin, sd=1/tau_emin)
  total_iterations$zVegMax[i] <- rtruncnorm(1, b=emax_high, mean=mu_emax, sd=1/tau_emax)
  total_iterations$zVegPeak[i]<- rtruncnorm(1, a=temp_emin, b=temp_emax, mean=mu_epeak, sd=1/tau_epeak)
  total_iterations$bMax[i] <- rtruncnorm(1, a=0, mean=mu_bmax, sd = 1/tau_bmax) / 10000 # g per m2 to g per cm2
  
  total_iterations$recalcitrantFrac[i] <-	rbeta(1, recalcitrant_alpha, recalcitrant_beta)
  total_iterations$rootDepthMax[i] <- rgamma(1, root_depth_alpha, reoot_depth_beta)
  total_iterations$omDecayRate[i] <- rbeta(1, decomp_alpha, decomp_beta)
  total_iterations$rootToShoot[i] <- rgamma(1, root_shoot_alpha, root_shoot_beta)
  total_iterations$rootTurnover[i] <- rgamma(1, root_turnover_alpha, root_turnover_beta)

  total_iterations$settlingVelocity[i] <- rgamma(1, settelingV_alpha, settelingV_beta)
  
  # Add to a table
  try({
  memCohortIteration <- runMemWithCohorts(startYear=2015, relSeaLevelRiseInit=0.3, rslrTotal=100,
                                        initElv=21.9, meanSeaLevel=7.4, meanHighWater=16.9, meanHighHighWater=25.4, meanHighHighWaterSpring=31.2, 
                                        suspendedSediment=3e-05, lunarNodalAmp=2.5, bMax=total_iterations$bMax[i], 
                                        zVegMin=total_iterations$zVegMin[i], zVegMax=total_iterations$zVegMax[i], zVegPeak=total_iterations$zVegPeak[i],
                                        plantElevationType="dimensionless", rootToShoot=total_iterations$rootToShoot[i],
                                        rootTurnover=total_iterations$rootTurnover[i], rootDepthMax=total_iterations$rootDepthMax[i], omDecayRate=total_iterations$omDecayRate[i],
                                        recalcitrantFrac=total_iterations$recalcitrantFrac[i], settlingVelocity=total_iterations$settlingVelocity[i],
                                        coreYear = 2050)
  
  # Calculate accumulation and sequestration
  # Save seq. rate
  total_iterations$OMAR[i] <- mean((memCohortIteration$core$slow_OM + 
                      memCohortIteration$core$fast_OM + 
                      memCohortIteration$core$root_mass) / memCohortIteration$core$input_yrs, 
                   na.rm=T) * 10000
  total_iterations$OMSR[i] <- mean(memCohortIteration$core$slow_OM / 
                               memCohortIteration$core$input_yrs, 
                             na.rm=T) * 10000
  })
  setTxtProgressBar(pb, i)
  
  
}

# Sensitivity analysis data frame
sensitivities <- data.frame(parameter=parameter_names,
                            type=parameter_types,
                          sensitivity=rep(NA, length(parameter_names)))

# Sensitivity analysis


for (j in 1:length(parameter_names)) {
  # For each parameter run on high setting
  for (k in 1:length(parameter_names)) {
    if (k == j) {
      assign(paste("temp_min_", parameter_names[k], sep=""), quantile(total_iterations[,parameter_names[k]],
                                                                      0.025,
                                                                      na.rm=T))
    } else {
      assign(paste("temp_min_", parameter_names[k], sep=""), quantile(total_iterations[,parameter_names[k]],
                                                                      0.5,
                                                                      na.rm=T))
    }
    if (k == j) {
      assign(paste("temp_max_", parameter_names[k], sep=""), quantile(total_iterations[,parameter_names[k]],
                                                                      0.975,
                                                                      na.rm=T))
    } else {
      assign(paste("temp_max_", parameter_names[k], sep=""), quantile(total_iterations[,parameter_names[k]],
                                                                      0.5,
                                                                      na.rm=T))
    }
  }
  # Run on high setting
  memCohortIteration_min <- runMemWithCohorts(startYear=2015, relSeaLevelRiseInit=0.3, rslrTotal=100,
                                          initElv=21.9, meanSeaLevel=7.4, meanHighWater=16.9, meanHighHighWater=25.4, meanHighHighWaterSpring=31.2, 
                                          suspendedSediment=3e-05, lunarNodalAmp=2.5, bMax=temp_min_bMax, 
                                          zVegMin=temp_min_zVegMin, zVegMax=temp_min_zVegMax, zVegPeak=temp_min_zVegPeak,
                                          plantElevationType="dimensionless", rootToShoot=temp_min_rootToShoot,
                                          rootTurnover=temp_min_rootTurnover, rootDepthMax=temp_min_rootDepthMax, omDecayRate=temp_min_omDecayRate,
                                          recalcitrantFrac=temp_min_recalcitrantFrac, settlingVelocity=temp_min_settlingVelocity,
                                          coreYear = 2050)
  seq_min <- mean(memCohortIteration_min$core$slow_OM / 
                                     memCohortIteration_min$core$input_yrs, 
                                   na.rm=T) * 10000
  # Run on low setting
  memCohortIteration_max <- runMemWithCohorts(startYear=2015, relSeaLevelRiseInit=0.3, rslrTotal=100,
                                              initElv=21.9, meanSeaLevel=7.4, meanHighWater=16.9, meanHighHighWater=25.4, meanHighHighWaterSpring=31.2, 
                                              suspendedSediment=3e-05, lunarNodalAmp=2.5, bMax=temp_max_bMax, 
                                              zVegMin=temp_max_zVegMin, zVegMax=temp_max_zVegMax, zVegPeak=temp_max_zVegPeak,
                                              plantElevationType="dimensionless", rootToShoot=temp_max_rootToShoot,
                                              rootTurnover=temp_max_rootTurnover, rootDepthMax=temp_max_rootDepthMax, omDecayRate=temp_max_omDecayRate,
                                              recalcitrantFrac=temp_max_recalcitrantFrac, settlingVelocity=temp_max_settlingVelocity,
                                              coreYear = 2050)
  seq_max <- mean(memCohortIteration_max$core$slow_OM / 
                    memCohortIteration_max$core$input_yrs, 
                  na.rm=T) * 10000
  # Calculate sensitivity
  # Store values
  sensitivities$sensitivity[j] <- abs(seq_max - seq_min)
  
}

sensitivities <- sensitivities %>%
  arrange(sensitivity)
sensitivities$order <- 1:nrow(sensitivities)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data=sensitivities, aes(x=order, y=sensitivity, color=type)) +
  geom_point(size=3, aes(shape=type, color=type)) +
  scale_x_discrete(limits = sensitivities$parameter) +
  xlab(element_blank()) +
  geom_segment(aes(x=order, xend=order, y=0, yend=sensitivity), size=1) +
  coord_flip() + 
  scale_colour_manual(values=cbbPalette) + 
  theme_minimal() + 
  labs(color=element_blank(),
       shape=element_blank())
ggsave("../../Desktop/MemSensitivitySeq.jpg", width=5, height=3.5)

