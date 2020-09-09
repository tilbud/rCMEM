#' Animate sediment cohorts accross an elevation transect
#' 
#' This function takes two runCohortMemTransect outputs, cohortsTransect and scenarioTransect tables, as inputs and visualizes soil formation as an animated .gif 
#' @param cohortsTransect data frame, annually tracked soil mass cohorts output from runCohortMemTransect
#' @param scenarioTransect data frame, annual summaries of inputs and outputs from runCohortMemTransect
#' @param filename, character, name of the output file
#' @param savePath character, filepath to save animation to
#' @param duration numeric, length in seconds of the animation
#' @param width numeric, width in inches of the .gif
#' @param height numeric, height in inches of the .gif
#' 
#' @export
animateCohortTransect <- function(scenarioTransect, cohortsTransect,
                                          duration = 30,
                                          width = 5,
                                          height = 3,
                                          savePath = getwd(),
                                          filename = "MEM-transect-animation.gif") {
  require(tidyverse, quietly = TRUE)
  require(gganimate, quietly = TRUE)
  require(gifski, quietly = TRUE)
  require(png, quietly = TRUE)
  
  surfaceElevations <- scenarioTransect %>% 
    dplyr::select(year, elvMax, elvMin, surfaceElevation) # %>% 
    # TODO rename years as year up stream so we don't have to rename it here
    # dplyr::rename(year = years)
  
  animateCohorts <- cohortsTransect %>% 
    dplyr::mutate(fraction_om = (fast_OM + slow_OM + root_mass) / 
                    (fast_OM + slow_OM + root_mass + mineral)) %>% 
    dplyr::select(year, layer_top, layer_bottom, elvMin, elvMax, fraction_om) %>% 
    dplyr::left_join(surfaceElevations, by = c("year", "elvMin", "elvMax")) %>% 
    dplyr::mutate(layer_top = surfaceElevation - layer_top,
                  layer_bottom = surfaceElevation - layer_bottom)
  
  waterLevel <- scenarioTransect %>% 
    dplyr::select(year, meanSeaLevel, meanHighWater) %>% 
    # dplyr::rename(year = years) %>% 
    group_by(year) %>%
    summarise(meanSeaLevel=first(meanSeaLevel), meanHighWater=first(meanHighWater)) %>% 
    gather(value="waterLevel", key="datum", -year) 
  
  soil_transect <- ggplot(data = animateCohorts, aes(xmin = elvMin, xmax = elvMax, ymin = layer_bottom, ymax = layer_top, 
                                                     frame = year
  )) +
    geom_rect(aes(fill = fraction_om), color = rgb(0,0,0, alpha = 0.1)) +
    scale_fill_gradient2(low = "darkgrey", mid = "lightgrey", high = "darkgreen", midpoint = 0.5, name = "Organic Matter (fraction)") + 
    theme_minimal() +
    geom_hline(data=waterLevel, aes(yintercept=waterLevel, lty=datum), color="blue") +
    ylab("Depth (cm NAVD88)") +
    xlab("Initial elevation (cm NAVD88)") +
    labs(title = 'Year: {round(frame_time,0)}') +
    transition_time(year) +
    ease_aes('linear')
  
  tempAnimation <- gganimate::animate(soil_transect, 
                                      duration = duration,
                                      nframes=length(unique(scenarioTransect$year)),
                                      renderer = gifski_renderer(),
                                      width = width, 
                                      height = height, 
                                      units = "in", 
                                      res = 300)
  (tempAnimation)
  
  gganimate::anim_save(filename=filename,
                       animation=tempAnimation,
                       path=savePath) 
}
