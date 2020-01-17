#' Depth of the given soil volume
#' 
#' Given a shape and mass for the roots, calculate the depth of a specific non-root volume. A specific solution is implimented for linear root volumes, and a general solution otherwise requires the \code{massLiveRoots.fn} to be specified. The general solution is a bineary search algorithm. 
#'
#' @param totalRootMass_per_area an integer that is the total mass per area of the roots, generally in g cm-3
#' @param rootDepthMax an integer that is the maximum root depth, generally in cm
#' @param rootDensity a numeric that is the root density in g per cm3
#' @param soilLength a numeric of the unit length of soil volume, generally 1 cm
#' @param soilWidth a numeric of the unitl witdth of the soil colume, generally 1cm
#' @param shape a character of the shape of the root distribution, currently only 'linear' is implimented
#' @param ... 
#' @param nonRootVolume.arr a numerical array of the cumulative non-rooting volume, assumed to be monotonically increasing
#' @param massLiveRoots.fn a function defining the mass of the live roots, if shape is not linear then this must be defined.
#' @param relTol a numeric the relative tolerence accepted for the general depth solution
#' @param verbose a boolean flag, true prints out a lot of strings to help debug the general solution
#'
#' @return A numeric corresponding to the depth of the specificed non root volume
#' 
#' @export
depthOfNonRootVolume <- function(nonRootVolume.arr, 
                                 massLiveRoots.fn = NULL,
                                 totalRootMass_per_area, 
                                 rootDepthMax, 
                                 rootDensity,
                                 shape = 'linear',
                                 soilLength = 1, soilWidth = 1,
                                 relTol = 1e-4,
                                 verbose = FALSE,
                                 ...){
  
  ####
  totalRootMass <- soilLength*soilWidth*totalRootMass_per_area

  if(verbose) print(paste('totalRootMass = ', totalRootMass))
  
  totalRootVolume <- totalRootMass/rootDensity
  if(verbose) print(paste('totalRootVolume = ', totalRootVolume))
  
  if(totalRootVolume > soilLength*soilWidth*rootDepthMax){
    stop('Bad root volume')
  }
  
  if(totalRootMass_per_area == 0){
    return(nonRootVolume.arr*soilLength*soilWidth)
  }
  
  ##If the shape is linear then solve the non-root volume depth explicitly 
  if(shape == 'linear'){
    if(verbose) print('We are linear.')
    rootWidth <- totalRootVolume*2/(rootDepthMax*soilLength)
    
    #nonRootVolume = ((rootWidth/rootDepthMax*depth^2)/2 + depth*(soilWidth-rootWidth))*soilLength
    #            0 = rootWidth / (2*rootDepthMax) * depth ^2 + 
    #                   (soilWidth-rootWidth) * depth - nonRootVolume/soilLength ##solve for depth
    coef1 <- rootWidth / (2*rootDepthMax)
    coef2 <- soilWidth-rootWidth
    coef3 <- -nonRootVolume.arr/soilLength
    ansDepth <- (-coef2 + sqrt(coef2^2-4*coef1*coef3))/(2*coef1)
    
    #correct for beyond root zone
    behondRootZone <- nonRootVolume.arr > soilLength*soilWidth*rootDepthMax - totalRootVolume
    
    if(verbose) print(paste0('nonRootVolume.arr: [', paste0(nonRootVolume.arr, collapse = ', '), ']'))
    if(verbose) print(paste0('nrv comparison: [', paste0(soilLength*soilWidth*rootDepthMax - totalRootVolume, collapse = ', '), ']'))
    if(verbose) print(paste0('beyondRootZon: [', paste0(behondRootZone, collapse = ', '), ']'))
    ansDepth[behondRootZone] <- (rootDepthMax +
                                   (nonRootVolume.arr - soilLength*soilWidth*rootDepthMax + totalRootVolume ) /
                                   (soilLength*soilWidth)) [behondRootZone] #treat as a square
    return(ansDepth)
    
  }else{ ##Otherwise apply a general binary search algorithm
    
    #What is the non rooting volumne down to rooting max
    nonRootVolumeToRootMax <- (soilLength*soilWidth*rootDepthMax) - totalRootVolume
    
    previousDepth <- 0
    possibleDepth.arr <- rep(NA, length=length(nonRootVolume.arr))
    
    for(nonRootVolume.index in (1:length(nonRootVolume.arr))){
      if(verbose)cat(paste('layerIndex: ', nonRootVolume.index, '----\n'))
      nonRootVolume <- nonRootVolume.arr[nonRootVolume.index]
      if(verbose)cat(paste('\ttarget nonRootVolume=', nonRootVolume, '\n'))
      if(nonRootVolume < relTol){
        if(verbose)cat('no volume\n')
        possibleDepth.arr[nonRootVolume] <- previousDepth
        next
      }
      #If we overfill that root zone volume
      if(nonRootVolumeToRootMax <= nonRootVolume){
        if(verbose)cat('outside root zone\n')
        #Subtract the root zone volumne, find the depth beyond and add to the max rooting depth
        possibleDepth.arr[nonRootVolume.index] <- (nonRootVolume-nonRootVolumeToRootMax)/(soilLength*soilWidth) + rootDepthMax
        previousDepth <- possibleDepth.arr[nonRootVolume.index]
        next
      }
      
      #incrament <- (min(rootDepthMax, nonRootVolume/(soilLength*soilWidth)) - previousDepth)
      #possibleDepth <- min(rootDepthMax, nonRootVolume/(soilLength*soilWidth)) - incrament / 2
      incrament <- rootDepthMax - previousDepth
      possibleDepth <- previousDepth + incrament/2
      
      maxSearchDepth <- ceiling(log(relTol, base=0.5)) #search down to relTol*depth
      for(ii in 1:maxSearchDepth){
        if(verbose)cat(paste('ii=', ii, '\n\tincrament=', incrament, '\n\tpossibleDepth=', possibleDepth))
        rootVolume <- massLiveRoots.fn(layerTop = 0, layerBottom = possibleDepth,
                                       totalRootMass_per_area = totalRootMass_per_area, 
                                       rootDepthMax=rootDepthMax, 
                                       soilLength=soilLength, soilWidth=soilWidth, shape=shape 
        ) / rootDensity
        if(verbose)cat(paste('\n\trootVolume_test=', rootVolume))
        nonRootVolume_test <- (possibleDepth * soilLength * soilWidth) - rootVolume
        if(verbose)cat(paste('\n\tnonRootVolume_test=', nonRootVolume_test))
        if( (abs(nonRootVolume-nonRootVolume_test) / nonRootVolume) < relTol){
          if(verbose)cat('\n\tdone\n')
          break
        }else{
          incrament <- incrament / 2
          
          if(verbose)cat(paste('\n\trelTol=', abs(nonRootVolume-nonRootVolume_test) / nonRootVolume))
          #Should we go up?
          if(nonRootVolume > nonRootVolume_test){
            if(verbose)cat('\n\tup\n')
            possibleDepth <- possibleDepth + incrament / 2
          }else{
            if(verbose)cat('\n\tdown\n')
            possibleDepth <- possibleDepth - incrament / 2
            
          }#if-else up/down
        }#if-else relTol
        
      }#for-loop refining search
      
      if(ii == maxSearchDepth){
        if(verbose) warning(paste('Target volumne precision [', relTol, '] is not reached'))
      }
      
      previousDepth <- possibleDepth
      possibleDepth.arr[nonRootVolume.index] <- possibleDepth
    }#for-loop across array
    
    return(possibleDepth.arr)
  }#if checking for unknown shape
  
}
