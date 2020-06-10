library(testthat)
context('precitedBiomass')
#testthat::test_file('tests/testthat/test-predictedBiomass.R')

test_that("predicted biomass gives default answers", {
  
  marshElevation <- 0
  bMax <- 2500
  zVegMax <- 3
  zVegMin <- -1
  
  agb <- predictedBiomass(0, 2500, 3, -1)
  
  zVegPeak <-(zVegMax+zVegMin)/2
  
  # From bmax, min, and max elevation limits, solve for parameters of a parabola.
  a <- -((-zVegMin * bMax - zVegMax * bMax) / ((zVegMin - zVegPeak) * (-zVegMax + zVegPeak)))
  b <- -(bMax / ((zVegMin - zVegPeak) * (-zVegMax + zVegPeak)))
  c <- (zVegMin * zVegMax * bMax) / ((zVegMin - zVegPeak) * (zVegMax - zVegPeak))
  
  
  expect_equal(agb, a*marshElevation + b*marshElevation^2 + c)
  #expect_equal(predictedBiomass(), agb)
})

test_that("check expected zero answers", {

  expect_equal(0, predictedBiomass(-4, 2500, 3, -1))
  expect_equal(0, predictedBiomass(5, 2500, 3, -1))
})

test_that('check stop conditions',{
  expect_error(predictedBiomass(0, 2500, -4, 1))
})
