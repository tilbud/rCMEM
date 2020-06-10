library(testthat)
context('depthOfNonRootVolumne')
#testthat::test_file('tests/testthat/test-depthOfNonRootVolumne.R')

test_that("test zero root mass", {
  
  expect_equal(depthOfNotRootVolume(nonRootVolume.arr = 1:4, 
                       massLiveRoots.fn = NULL,
                       totalRootMassPerArea = 0, 
                       rootDepthMax = 2, 
                       rootDensity = 1.3,
                       shape = 'linear',
                       soilLength = 1, soilWidth = 1,
                       relTol = 1e-4,
                       verbose = FALSE), 1:4)
})