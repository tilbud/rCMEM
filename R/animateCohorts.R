#' Animate sediment cohorts
#' 
#' This function takes two MEM outputs, cohorts and scenario tables, as inputs and visualizes soil formation as an animated .gif 
#' @param cohorts data frame, annually tracked soil mass cohorts output from runMemWithCohorts
#' @param scenario data frame, annual summaries of inputs and outputs from runMemWithCohorts
#' @param filename, character, name of the output file
#' @param savePath character, filepath to save animation to
#' @param chPallette vector, a vector of colors to use to symbolize the different mass cohorts
#' @param trackThresholds vector, a vector of characters indicating which water leves in the scenario table to map as horizontal lines
#' @param duration numeric, length in seconds of the animation
#' @param width numeric, width in inches of the .gif
#' @param height numeric, height in inches of the .gif
#' 
#' @export
animateCohorts <- function(cohorts, scenario,
                           filename = "MEM-CTM-animated.gif",
                           savePath = getwd(),
                           chPalette = c("#56B4E9", "#999999", "#E69F00", "#009E73"), 
                           trackThresholds = c("MSL", "MHW"), duration = 30,
                           width = 4.5, height = 4.5) {
  
  require(tidyverse, quietly = TRUE)
  require(gganimate, quietly = TRUE)
  require(gifski, quietly = TRUE)
  require(png, quietly = TRUE)
  
  surface_elv <- scenario %>%
    dplyr::select(years, surfaceElevation) %>%
    dplyr::rename(year=years)
  
  # First reshape the mass cohorts so that they're in long form
  mass_cohorts <- cohorts %>%
    dplyr::select(-cumCohortVol, -respired_OM) %>% 
    tidyr::gather(key = "mass_pool", value = "mass_fraction", 
           -age, -year, -layer_top, -layer_bottom) %>%
    dplyr::group_by(year, age) %>%
    dplyr::mutate(mass_pool = factor(mass_pool, 
                              levels=c("mineral", 
                                       "slow_OM", 
                                       "fast_OM", 
                                       "root_mass"))) %>%
    dplyr::arrange(year, age, mass_pool) %>%
    dplyr::mutate(max_mass = cumsum(mass_fraction),
           min_mass = ifelse(mass_pool==first(mass_pool),0,lag(max_mass)),
           mass_pool = as.character(mass_pool)) %>%
    # Join mass cohorts with scenario table to convert depths to referenced elevations
    dplyr::ungroup() %>%
    dplyr::left_join(surface_elv, by="year") %>%
    dplyr::mutate(layer_top=surfaceElevation-layer_top, 
           layer_bottom=surfaceElevation-layer_bottom)
  
  # Reshape the scenario table
  tides <- scenario %>%
    # Track any elevation threholds in the animation speciefied.
    # MSL and MHW are the defaults
    dplyr::select(years, trackThresholds) %>%
    tidyr::gather(value="WaterLevel", key="datum", -years) %>%
    dplyr::rename(year=years) %>%
    dplyr::arrange(year) %>%
    dplyr::filter(complete.cases(.))
  
  # get rid of any NA values.               
  mass_cohorts_almostAll <- mass_cohorts %>%  dplyr::filter(complete.cases(.))
  
  # gganimate stuff
  animate_mass_cohorts <- ggplot2::ggplot(data = mass_cohorts_almostAll, 
                                 aes(xmin = min_mass, xmax = max_mass, 
                                     ymin = layer_top, ymax = layer_bottom , 
                                     frame = year
                                 )) +
    ggplot2::geom_rect(aes(fill = mass_pool), color = rgb(0,0,0, alpha = 0.1)) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values=chPalette) +
    ggplot2::geom_hline(data=tides, aes(yintercept=WaterLevel, lty=datum), color="blue") +
    ggplot2::theme(text = element_text(size=14)) +
    ggplot2::ylab("Depth (cm)") +
    ggplot2::xlab("Mass Accumulated Per Cohort (g)") +
    ggplot2::labs(title = 'Year: {round(frame_time,0)}',
         fill = "Mass Pools") +
    gganimate::transition_time(year) +
    gganimate::ease_aes('linear')
  
  tempAnimation <- gganimate::animate(animate_mass_cohorts, 
                                      duration = duration,
                                      renderer = gifski_renderer(),
                                      width = width, 
                                      height = height, 
                                      units = "in", 
                                      res = 300)
  (tempAnimation)
  # save gif to filepath
  gganimate::anim_save(filename=filename,
                       animation=tempAnimation,
                       path=savePath)
}
