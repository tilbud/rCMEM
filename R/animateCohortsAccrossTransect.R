
animateCohortsAccrossTransect <- function(scenarioTransect, cohortsTransect,
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
    dplyr::select(years, elvMax, elvMin, surfaceElevation) %>% 
    dplyr::rename(year = years)
  
  animateCohorts <- cohortsTransect %>% 
    dplyr::mutate(fraction_om = (fast_OM + slow_OM + root_mass) / 
                    (fast_OM + slow_OM + root_mass + mineral)) %>% 
    dplyr::select(year, layer_top, layer_bottom, elvMin, elvMax, fraction_om) %>% 
    dplyr::left_join(surfaceElevations, by = c("year", "elvMin", "elvMax")) %>% 
    dplyr::mutate(layer_top = surfaceElevation - layer_top,
                  layer_bottom = surfaceElevation - layer_bottom)
  
  waterLevel <- scenarioTransect %>% 
    dplyr::select(years, MSL, MHW) %>% 
    dplyr::rename(year = years) %>% 
    group_by(year) %>%
    summarise(MSL=first(MSL), MHW=first(MHW)) %>% 
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
