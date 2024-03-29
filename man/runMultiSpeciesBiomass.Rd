% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runMultiSpeciesBiomass.R
\name{runMultiSpeciesBiomass}
\alias{runMultiSpeciesBiomass}
\title{Run Multiple Species Biomass Functions}
\usage{
runMultiSpeciesBiomass(
  z,
  bMax,
  zVegMax,
  zVegMin,
  zVegPeak = NA,
  rootToShoot,
  rootTurnover,
  abovegroundTurnover = NA,
  rootDepthMax,
  speciesCode = NA,
  competition.fn = NA
)
}
\arguments{
\item{z}{a numeric, elevation}

\item{bMax}{a numeric, or vector of numerics, maximum biomass}

\item{zVegMax}{a numeric, or vector of numerics, upper elevation of biomass limit}

\item{zVegMin}{a numeric, or vector of numerics, lower elevation of biomass limit}

\item{zVegPeak}{(optional) a numeric, or vector of numerics, elevation of peak biomass}

\item{rootToShoot}{a numeric, or vector of numerics, root to shoot ratio}

\item{rootTurnover}{a numeric, or vector of numerics, belowground biomass annual turnover rate}

\item{abovegroundTurnover}{(optional) a numeric, or vector of numerics, aboveground biomass annual turnover rate}

\item{rootDepthMax}{a numeric, or vector of numerics, maximum (95\%) rooting depth}

\item{speciesCode}{(optional) a character, or vector of characters, species names or codes associated with biological inputs}

\item{competition.fn}{(optional) a function that takes at least a dataframe bio_table as an input, models competition between multiple species, and outputs an one row data frame aggregating biomass and biological inputs}
}
\value{
a one dataframe with above ground biomass, and biological parameters representing the dominant specie(s) at the elevation
}
\description{
This function takes elevation and biological inputs, runs parabolic predictBiomass function for multiple species, and returns a single above ground biomass value and set of biological parameters based on a competition function.
}
