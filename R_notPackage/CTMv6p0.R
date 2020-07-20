# This file will house the functions for translating CTM v 6.1 into R functions

# Inputs and key parameters
{
  Dm <-	30 # Depth of root zone (cm)
  RT <-	1000 # Root rhizome biomass (g/m^2)
  kr <-	0.2 # Refractory fraction of OM production
  SSC <-	10 # Suspended Sediment Concentration (mg/L)
  OMDR <-	0.8 # Labile OM decay rate (yr^-1)
  BGTR <-	0.5 # Below ground turnover rate (yr^-1) 
  D <-	5 # Depth of surface below meanHighWater (cm^3)
  CY <-	2015 # Core year
  k1 <-	0.085 # OM self packing density (g/cm^3)
  k2 <-	1.99 # Mineral sefl packing density (g/cm^3)
  expRoot <- F
}

# Create vectors for values we will store for subannual cohorts
{
  R <- c() # Root mass
  Rvol <- c() # Root volume
  mineral <- c() # mineral mass
  LabileOM <- c() # Labile Organic matter
  refracOM <- c() # Refractory Organic matter
  yrs <- c() # Years
  cumage <- c() # cumulative age
  wtage <- c() # weighted age
  carbagewt <- c() # carbon weighted age
  LOI <- c() # loss on ignition
  yrsold <- c() # years old
  decay <- c() # decay
  Labprod <- c() # labile production
  labsec <- c() # Labile section
  Pb210 <- c() # 210 Pb
  totvol <- c() # total volume
  Minvol <- c() # mineral volume
  Rootvol <- c() # root volume
  Slowvol <- c() # Slow decay om pool volume
  Fastvol <- c() # Fast decay om pool volume
  cumvol <- c() # Cumulative volume
  iterationsToSolve <- c()
}

# Make conversions and constants
{
  DT <- 0.1 # Time step (1 / 10th of a year)
  m <- SSC * 0.000001
  d <- D
  n_tides <- 704
}

# Convert rate constants to subannual scale
{
  # These rate constants are multiplied by DT (DT is 1/10 of a year)
  vol1 <- BGTR * kr * 0.0001 * RT / k1 + d * m * n_tides / k2 
  
  # If the steady state volume production is greater than 1 cm3
  # ... then lower the time step to 1/100th of a year.
  if (vol1 > 1) { DT <- 0.01 }
  
  # Convert BGTR from annual to subtep of the algorithm
  bgtrSubAnn <- BGTR * DT
  
  # m is technically not a rate constant, 
  # but it is used in the formula for rate of mineral deposition
  mSubAnn <- m * DT
  
  # multiply the decay rate by subannual time step
  omdrSubAnn <- OMDR * DT
}

# Set Cs peak year
{
  # 1963 is the year of the Cs peak, its position or age in the simulation
  # depends on the starting year CY
  Csyr <- CY - 1963
  perm <- F # Set
}

# Define steady state accretion rate
Accrate <- (kr * BGTR * RT * 0.0001) / k1 + 704 * d * m / k2

# Begin the interation over 10000 years

# Add the refracOM from surface cohorts
# Compute the decay fraction of labile OM in cohort i
# = sum of (1-omdr) + (1-omdr)^2 + (1-omdr)^n where n is the cohort number
# cycle through y2 to find y2=y1

for (i in 1:10000) {
  # Tthe mineral volume is constant over time
  Minvol <- c(Minvol, n_tides * mSubAnn * d / k2) # !!! When I come back to this replace the loop with a rep
}

# THIS SECTION IS FOR THE EXPONENTIAL ROOT MODEL

{
  if (expRoot == T) { # !!! Replace this conditional with function based around root shape
  
    # !!! replace solver with a function
    # !!! better comment how solver works.
    # !!! Make this first step of a function
    # !!! Conditionality, if this is first first iteration / surface
    
    kd <- log(0.05) / Dm
    Rm <- -0.95 * RT * 0.0001 * kd / (1 - exp(kd * Dm))
    
    d2 <- exp(kd * Minvol[1]) # d2 is the greater depth
    
    Rvol[1] <- (Rm / kd) * (d2 - 1) / k1 # live roots in first mineral section
    Fastvol[1] <- Rvol[1] * bgtrSubAnn * (1 - kr) * (1 - omdrSubAnn) # labile OM volume
    Slowvol[1] <- Rvol[1] * bgtrSubAnn * kr # refractory OM volume
    totvol[1] <- Minvol[1] + Rvol[1] + Fastvol[1] + Slowvol[1]
    cumvol[1] <- totvol[1]
    
    d2_temp <- c(d2)
    Rvol_temp <-c(Rvol[1])
    totvol_temp <-c(totvol[1])
    d2change_temp <-c(NA)
    
    # This section converges on a solution
    k <- 1
    stop <- F
    while ((stop == F) & (k < 20)) {
      
      # !!! replace solver with a function
      # !!! better comment how solver works.
      # !!! Make this second step of a function
      # !!! Conditionality, if this is first iteration / surface
      
      # between 0 and Dm = (Ro/c)[exp(cDm) - exp(c0)]
      d2 <- exp(kd * cumvol[1]) # d2 is the greater depth
      Rvol[1] <- (Rm / kd) * (d2 - 1) / k1 # live root vol in first mineral section
      Slowvol[1] <- Rvol[1] * bgtrSubAnn * kr
      Fastvol[1] <- Rvol[1] * bgtrSubAnn * (1 - kr) * (1 - omdrSubAnn)
      totvol[1] <- Minvol[1] + Rvol[1] + Fastvol[1] + Slowvol[1]
      cumvol[1] <- totvol[1]
      d2new <- exp(kd * cumvol[1])
      if (abs(d2 - d2new) < 0.0000001) { stop <- T }
      J <- k
      k <- k + 1
      
      d2_temp <- c(d2_temp, d2)
      Rvol_temp <-c(Rvol_temp, Rvol[1])
      totvol_temp <-c(totvol_temp, totvol[1])
      d2change_temp <-c(d2change_temp, abs(d2 - d2new))
      
    }
    
    for (i in 2:5999) {
      
      # !!! replace solver with a function
      # !!! better comment how solver works.
      # !!! Make this third step of a function
      # !!! Conditionality, if this is not first iteration / surface
      
      d1 <- exp(kd * (cumvol[i - 1])) # the top of the section, closer to the surface
      d2 <- exp(kd * (cumvol[i - 1] + totvol[i - 1])) # d2 is a greater depth
      Rvol[i] <- (Rm / kd) * (d2 - d1) / k1 # live roots in first
      
      # mineral section first guess
      Fastvol[i] <- (Rvol[i] * bgtrSubAnn * (1 - kr) + Fastvol[i - 1]) * (1 - omdrSubAnn)
      Slowvol[i] <- (Rvol[i] * bgtrSubAnn * kr + Slowvol[i - 1])
      totvol[i] <- Minvol[i] + Rvol[i] + Fastvol[i] + Slowvol[i]
      cumvol[i] <- cumvol[i - 1] + totvol[i]
      
      d2_temp <- c(d2)
      Rvol_temp <-c(Rvol[i])
      totvol_temp <-c(totvol[i])
      d2change_temp <-c(NA)
      
      # This section converges on a solution
      k <- 1
      stop <- F
      while ((stop == F) & (k < 20)) {
        
        # !!! replace solver with a function
        # !!! better comment how solver works.
        # !!! Make this fourth step of a function
        # !!! Conditionality, if this is not first iteration / surface
        
        d2 <- exp(kd * (cumvol[i - 1] + totvol[i])) # adjust volume to accomdate roots
        Rvol[i] <- (Rm / kd) * (d2 - d1) / k1 # live roots in first
        
        # mineral section
        Slowvol[i] <- Rvol[i] * bgtrSubAnn * kr + Slowvol[i - 1]
        Fastvol[i] <- (Rvol[i] * bgtrSubAnn * (1 - kr) + Fastvol[i - 1]) * (1 - omdrSubAnn)
        totvol[i] <- Minvol[i] + Rvol[i] + Fastvol[i] + Slowvol[i]
        d2new <- exp(kd * (cumvol[i - 1] + totvol[i]))
        
        if (abs(d2 - d2new) < 0.0000001) { stop <- T }
       
        J <- k
        k <- k + 1
        
        cumvol[i] <- cumvol[i - 1] + totvol[i]
        
        d2_temp <- c(d2_temp, d2)
        Rvol_temp <-c(Rvol_temp, Rvol[i])
        totvol_temp <-c(totvol_temp, totvol[i])
        d2change_temp <-c(d2change_temp, abs(d2 - d2new))
      }
    }
  } else { # If the root profile is linear
      
      # !!! Replace this conditionality with a root input function
      
      # !!! A lot of this will be simplified with a single solver function
      
      d2 <- Minvol[1]
      Rm <- 0.0001 * RT / (0.5 * Dm)
      kd <- RT * 0.0001 / (0.5 * Dm ^ 2)
      Rvol[1] <- (Rm * d2 - (kd / 2) * (d2 ^ 2)) / k1
      Fastvol[1] <- Rvol[1] * bgtrSubAnn * (1 - kr) * (1 - omdrSubAnn) # labile OM volume
      Slowvol[1] <- Rvol[1] * bgtrSubAnn * kr # refractory OM volume
      totvol[1] <- Minvol[1] + Rvol[1] + Fastvol[1] + Slowvol[1]
      cumvol[1] <- totvol[1]
      
      d2_temp <- c(d2)
      Rvol_temp <-c(Rvol[1])
      totvol_temp <-c(totvol[1])
      d2change_temp <-c(NA)
      
      # This section converges on a solution
      k <- 1
      stop <- F
      while ((stop == F) & (k < 20)) {
        
        # !!! A lot of this will be simplified with a single solver function
        
        d2 <- cumvol[1] # adjust volume to accomdate roots
        Rvol[1] <- (Rm * d2 - (kd / 2) * (d2 ^ 2)) / k1 # live roots in first
        
        # mineral section
        Slowvol[1] <- Rvol[1] * bgtrSubAnn * kr
        Fastvol[1] <- Rvol[1] * bgtrSubAnn * (1 - kr) * (1 - omdrSubAnn)
        totvol[1] <- Minvol[1] + Rvol[1] + Fastvol[1] + Slowvol[1]
        cumvol[1] <- totvol[1]
        d2new <- cumvol[1]
        
        if (abs(d2 - d2new) < 0.0000001) { stop <- T }
        
        J <- k # !!! Not really sure what J is
        k <- k + 1
        
        d2_temp <- c(d2_temp, d2)
        Rvol_temp <-c(Rvol_temp, Rvol[1])
        totvol_temp <-c(totvol_temp, totvol[1])
        d2change_temp <-c(d2change_temp, abs(d2 - d2new))
        
      }
      
      # !!! A lot of this seems like it was copied from an old sheet because it's out of place here
      
      # refracOM[1] <- Slowvol[1] * k1
      # sumrootprod <- Rvol[1] * bgtrSubAnn
      # decay[1] <- (Rvol[1] * bgtrSubAnn * (1 - kr)) * omdrSubAnn # annual decay
      # sumanndecay <- sumanndecay + decay[1]
      # LabileOM[1] <- Fastvol[1] * k1
      # Pb210[1] <- Pb210surf * Exp(-DT * 0.0311)
      # volcheck <- volcheck + (refracOM[1] + R[1] + LabileOM(1)) / k1 +
      #    mineral(1) / k2
      # LOI[1] = 90 * (refracOM[1] + LabileOM[1] + R[1]) / (refracOM[1] + LabileOM[1] + R[1] + mineral[1])
      
      for (i in 2:5999) {
        d1 <- cumvol[i - 1]
        d2 <- d1 + totvol[i - 1] # to initialize, use the last total volume
        test <- d2 - d1 # !!!  don't think I need this
        Rvol[i] <- max(0, (Rm * d2 - (kd / 2) * (d2 ^ 2) - ((Rm * d1 - (kd / 2) * (d1 ^ 2)))) / k1)
        Fastvol[i] <- (Rvol[i] * bgtrSubAnn * (1 - kr) + Fastvol[i - 1]) * (1 - omdrSubAnn)
        Slowvol[i] <- (Rvol[i] * bgtrSubAnn * kr + Slowvol[i - 1])
        totvol[i] <- Minvol[i] + Rvol[i] + Fastvol[i] + Slowvol[i]
        cumvol[i] <- cumvol[i - 1] + totvol[i]
        
        # This section converges on a solution
        k <- 1
        stop <- F
        while ((stop == F) & (k < 20)) {
          
          d2 <- cumvol[i - 1] + totvol[i]
          Rvol[i] <- max(0, (Rm * d2 - (kd / 2) * (d2 ^ 2) - ((Rm * d1 - (kd / 2) * (d1 ^ 2)))) / k1)
          Slowvol[i] <- Rvol[i] * bgtrSubAnn * kr + Slowvol[i - 1]
          Fastvol[i] <- (Rvol[i] * bgtrSubAnn * (1 - kr) + Fastvol[i - 1]) * (1 - omdrSubAnn)
          totvol[i] <- Minvol[i] + Rvol[i] + Fastvol[i] + Slowvol[i]
          d2new <- cumvol[i - 1] + totvol[i]
          
          if (abs(d2 - d2new) < 0.0000001) { stop <- T }
          
          cumvol[i] <- cumvol[i - 1] + totvol[i]
          
          J <- k
          k <- k + 1
          
          d2_temp <- c(d2_temp, d2)
          Rvol_temp <-c(Rvol_temp, Rvol[i])
          totvol_temp <-c(totvol_temp, totvol[i])
          d2change_temp <-c(d2change_temp, abs(d2 - d2new))
          
        }
      }
    }
}

# SUMMARIZE PROFILES
#  W.J.M. van der Linden, S.A.P.L. Cloetingh, ?J.P.H. Kaasschieter - 2013 Science
# The average amount of Pb210 transferred from the water onto the suspended sediment in the Waddensea was found to be
# 0.001397 dpm/cm2/day 
# (mean input from the atmosphere being 0.0014 dpm/cm2/day = .51 dpm/cm2/yr)

{
  Pb210surf <- 0.51 # dpm cm^-2 yr^-1 with 22.3 yr halflife
  sumrootprod <- 0 # annual root production
  refracprod <- 0 # annual refractor production
  sumanndecay <- 0
  totroot <- 0
  minsum <- 0
  labsum <- 0
  refsum <- 0
  rootsum <- 0
  refracprod <- 0
  decaysum <- 0
  Pbsum <- 0
  totroot <- 0
  volcheck <- 0
  Line <- 3
  iend <- 0
  section <- 1
}

# !!! can probably replace this whole loop with dplyr operations
# !!! may be faster to do things in base R although less

for (i in 1:5999) {
  R[i] <- Rvol[i] * k1 # total root biomass in annual section
  refracOM[i] <- Slowvol[i] * k1 # total refraction production in section
  Labprod[i] <- Fastvol[i] * k1 # annual production of labile in section
  LabileOM[i] <- Fastvol[i] * k1 # !!! redundent.
  mineral[i] <- Minvol[i] * k2 # g/cm2
  volcheck <- volcheck + (refracOM[i] + R[i] + LabileOM[i]) / k1 + mineral[i] / k2
  cumage[i] <- i
  
  # analytical decay = (1 - kr) * BGTR * RT
  # !!! this calculation doesn't work in R yet for i = 1
  decay[i] <- (R[i] * bgtrSubAnn * (1 - kr) + LabileOM[i - 1]) * (omdrSubAnn) # annual decay
  
  LOI[i] <- 90 * (refracOM[i] + LabileOM[i] + R[i]) / 
    (refracOM[i] + LabileOM[i] + R[i] + mineral[i])
  
  Total <- refracOM[i] + LabileOM[i]
  Pb210[i] <- Pb210surf * exp(-DT * 0.0311 * i)
  
  # analytical root prod is RT*kr*BGTR
  totroot <- totroot + R[i]
  sumrootprod <- sumrootprod + R[i] * bgtrSubAnn # annual root production
  refracprod <- refracprod + kr * bgtrSubAnn * R[i] # annual refractor production
  sumanndecay <- sumanndecay + decay[i] # annual decay
  Pbsum <- Pbsum + Pb210[i] # total counts (dpm) per cm slice
  minsum <- minsum + mineral[i]
  labsum <- labsum + LabileOM[i]
  refsum <- refsum + refracOM[i]
  rootsum <- rootsum + R[i]
  decaysum <- decaysum + decay[i]
  
  if (i * DT == Csyr) {
    secCsyr <- i * DT
    Csdepth <- cumvol[i]
    acrCsyr <- cumvol[i] / secCsyr
  }
  
  # End If from test of section is equal to 1 cc of volume
  
  if (cumvol[i] >= section) {
    
    # !!! need to figure out the functionality of this line and document it better
    yrs[section] <- (i - iend) * DT # calc the number of years in a section
    
    
    if (iend == 0) { # !!! what is iend?
      
      yrszero <- yrs[section] 
      
      # outputs.Cells(Line, 1) = cumvol[i]
      # outputs.Cells(Line, 2) = yrs(section)
      
      #  If exproot Then
      #  If iend > 0 Then Rbiom = 10000 * (Rm / kd) * (Exp(kd * cumvol[i]) -
      #  Exp(kd * cumvol(iend))) ' root biomass
      #  If iend = 0 Then Rbiom = 10000 * (Rm / kd) * (Exp(kd * cumvol[i]) -
      #  1) ' root biomass
    
      } else {

      #  Else:
      
      #  If iend > 0 Then Rbiom = 10000 * WorksheetFunction.Max(0, (Rm *
      #  cumvol[i] - (kd / 2) * (cumvol[i] ^ 2) - (Rm * cumvol(iend) - (kd / 2) *
      #    (cumvol(iend) ^ 2)))) ' root production per year
      
      #  If iend = 0 Then Rbiom = 10000 * WorksheetFunction.Max(0, (Rm - Rm *
      #    (cumvol[i] ^ 2) / (2 * Dm))) ' root production per year
      
      #  End If
    }
      
      #     outputs.Cells(Line, 10) = BGTR * Rbiom ' root production per year
      #     outputs.Cells(Line, 10) = BGTR * Rbiom ' root production per year
      #                                            
      #     outputs.Cells(Line, 5) = Rbiom ' root biomass
      #     outputs.Cells(Line, 6) = 10000 * labsum
      #     outputs.Cells(Line, 7) = 10000 * refsum
      #     outputs.Cells(Line, 8) = 10000 * minsum ' mineral
      #     drywt = 0.0001 * Rbiom + labsum + refsum + minsum
      #     LossOnIgnition = 100 * (0.0001 * Rbiom + labsum + refsum) / drywt
      #     outputs.Cells(Line, 12) = LossOnIgnition ' LOI
      #                                            
      #     outputs.Cells(Line, 13) = 10 * (cumvol[i] - cumvol(iend)) / yrs(section)
      #     ' vertical accretion per section mm/yr
      #     'outputs.Cells(line, 13) = 10 * outputs.Cells(line, 1) /
      #     outputs.Cells(line, 2) ' accretion per section
      #     outputs.Cells(Line, 3) = volcheck
      #     outputs.Cells(Line, 4) = I * DT 'cumage[i]
      #     outputs.Cells(Line, 9) = 10000 * decaysum / yrs(section)
      #                                            
      #     'Pb210[i] = Pb210surf * Exp(-i * DT * 0.69315) / minsum
      #                                            
      #     Pb210(section) = Pb210surf * yrszero * Exp(-DT * I * 0.5 * 0.69315) / minsum
      #                                            
      #                                            
      #     outputs.Cells(Line, 11) = 1 / (0.01 * LossOnIgnition / k1 + (1 - 0.01 *
      #     LossOnIgnition) / k2) ' bulk density
      #     'outputs.Cells(line, 11) = 1 / (0.01 * LOI(section) / k1 + (1 - 0.01 *
      #     LOI(section)) / k2)
      #     outputs.Cells(Line, 15) = Pbsum / minsum
      #     outputs.Cells(Line, 14) = Log(Pbsum / minsum)
      #                                            
      #     Line = Line + 1
      #     section = cumvol[i] + 1
      #     iend = I
      #     minsum = 0
      #     labsum = 0
      #     refsum = 0
      #     rootsum = 0
      #     labsum = 0
      #     decaysum = 0
      #     volcheck = 0
      #     Pbsum = 0
      #     volcheck = volcheck + (refracOM(1) + R(1) + LabileOM(1)) / k1 +
      #     mineral(1) / k2
      #     LOI(1) = 90 * (refracOM(1) + LabileOM(1) + R(1)) / (refracOM(1) +
      #     LabileOM(1) + R(1) + mineral(1))
      #     '**********************************************
      #     End If
    }
    # Next i
  }
  
  # 'The first model, known as CF:CS (the Constant Flux Constant
  # 'Sedimentation Rate model), assumes that there is a constant flux of 210Pb
  # 'and that the rate of sediment deposition is constant as well. With this
  # model
  # 'the sedimentation rate can be calculated using the slope of the line
  # derived
  # 'from the linear regression of ln210Pbex and the depth layer according
  # to the
  # 'following equations (Bierman et al. 1998):
  #   'Ax = A0exp(-bx) , v =lamda/b, where Ax – excess 210Pb activity at depth
  # x [Bq kg-1 d.m.], A0 is activity at
  # 'the surface layer [Bq kg-1 d.m.], b is the slope defined by regression
  # through
  # 'the data, x is depth [cm], v is sedimentation rate – LAR [cm year-1]
  # and lamda is
  # 'the 210Pb radioactive decay constant (0.03114 year-1).
  # ' ln(Ax)=ln(Ao)-bx
  # ' The CRS model was initiated by Appleby and
  # 'Oldfield (1978) and Oldfield and Appleby (1984) assumes a constant
  # 210Pb flux at water-sediment interface
  # 'and requires both the integrated activity and the differential activity
  # to yield a variable sedimentation rate
  # '(Shukla 2002). The CRS model is used when the supply of 210Pbex is
  # constant and the sediment deposition
  # 'rate is variable (Appleby and Oldfield 1978; Robbins 1978). The CIC
  # model has been applied by Robbins
  # 'and Edgington (1975), Durham and Joshi (1980) among others assumes
  # both a constant mass flux and a
  # 'constant activity flux at the water-sediment interface; and
  # requires only differential activity to yield a
  # 'constant sedimentation rate (Shukla 2002). It also appropriate when
  # initial activity of 210Pbex is constant and
  # 'there is no mixing of surface sediments (Robbins and Edgington
  #                                                  1975). This implies deposition of
  # 'sedimentary material characterized by constant 210Pbex activity such
  # that either both
  # '210Pbex activity and mass of deposited material to the sediment
  # surface are constant or both vary at the same
  # 'rate (Zaborska et al. 2007).
  # 'Pb210(1) = Pb210surf * Exp(-1 * 0.69315 / 22.3) / mineral(1)
  # 
  # 'Arg1 = outputs.Cells.Range("d3:d82")
  # 
  # 'Arg2 = outputs.Cells.Range("o3:o82")
  # Pbslope = instance.Slope(outputs.Cells.Range("n3:n46"),
  #                          outputs.Cells.Range("a3:a46"))
  # Dhalf = -0.69315 / Pbslope
  # '=-0.03114/((LN(G30/G2))/(B30-B2))*10 Craft formula
  # PBAR = Dhalf / 22.3
  # '*************************************** end of loop
  # ******************************************
}

# surfacr = totvol(1) / DT ' cm
# Cells(3, 12) = "Steady St Accretion"
# Cells(3, 13) = 10 * (kr * BGTR * RT * 0.0001 / k1 + 704 * d * m / k2) /
# DT ' mm per year accretion
# Cells(3, 14) = "mm/yr"
# Cells(4, 12) = "Surface Accretion"
# Cells(5, 12) = "Accretion to Cs pk"
# Cells(7, 12) = "Weighted C Age"
# Cells(8, 12) = "Cesium Horizon Age"
# Cells(8, 14) = " @ " & Round(Csdepth, 1) & " cm"
# Cells(4, 13) = surfacr * 10
# Cells(5, 13) = acrCsyr * 10
# Cells(6, 13) = PBAR * 10
# Cells(7, 14) = "yr"
# ' numerical calculations
# Cells(4, 7) = sumanndecay * 10000 / DT ' annual production
# Cells(4, 8) = sumrootprod * 10000 / DT ' annual production
# Cells(4, 9) = Cells(4, 8) * kr ' annual production of refractory
# material (g/m2/yr)
# ' analytical calculations
# Cells(5, 7) = (1 - kr) * BGTR * RT / DT ' annual decay rate (g/m2/yr)
# Cells(5, 8) = BGTR * RT / DT ' annual production
# Cells(5, 9) = kr * BGTR * RT / DT
# 
# ' The organic carbon at the depth of the Cs peak is a variety of ages.
# ' This section computes the weighted age of carbon at the peak Cs horizon.
# For I = 1 To secCsyr
# yrsold[i] = 0
# ' yrs(j) is the number of years of inputs in each 1 cm section
# For J = I To secCsyr
# yrsold[i] = yrsold[i] + yrs(J)
# Next J
# Next I
