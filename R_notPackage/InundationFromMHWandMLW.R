
calcTimeInundated <- function(elevation, meanHighWater=0.364, meanLowWater=0.067) {
  # Equation from in review paper - Tidal inundation modeling within GIS
  # Sent to JH by Bob Hickey 12 March 2019.
  # (Adapted from Department of the Navy. Hydrographic Branch, 1994)
  if (elevation > meanHighWater) {
    inundationTime <- 0
  } else if (elevation < meanLowWater) {
    inundationTime <- 6.21*2
  } else {
    # Rising time over cell = 6.21 (A/pi - 1)
    # A = 2*pi - cos-1 [2 (height of cell – meanLowWater) / (meanHighWater – meanLowWater) - 1] radians
    A1 <- 2*pi-acos(2 * (elevation - meanLowWater) / (meanHighWater-meanLowWater) - 1)
    risingTime <- 6.21 * (A1 / pi - 1)
    # A = 2*pi - cos-1 [2 (height of cell – meanHighWater) / (meanLowWater – meanHighWater) - 1] radians
    A2 <- 2*pi-acos(2 * (elevation - meanHighWater) / (meanLowWater-meanHighWater) - 1)
    # Falling time over cell = 6.21 (A/p - 1) where
    fallingTime <- 6.21 * (A2 / pi - 1)
    # Cell inundation time = abs (time rising - 6.21) + time falling
    inundationTime <- abs(risingTime - 6.21) + fallingTime
  }
} 

tideElevationRange <- seq(meanLowWater-0.03, meanHighWater+0.03, by = 0.005)

oneTidalCycle <- sapply(tideElevationRange, calcTimeInundated)

plot(tideElevationRange, oneTidalCycle, type="l")
