library(tidyverse)

download.file(url = "https://tidesandcurrents.noaa.gov/api/datagetter?product=monthly_mean&application=NOS.COOPS.TAC.WL&station=8575512&begin_date=19280901&end_date=20180101&datum=NAVD&units=metric&time_zone=GMT&format=csv",
              destfile = "vignettes/sampleData/annapolisTideGauge.csv")

annapolisTideGauge <- read_csv("vignettes/sampleData/annapolisTideGauge.csv")

annapolisTideGaugeTimeSeries <- annapolisTideGauge %>%
  group_by(Year) %>%
  summarise(MHW = mean(MHW, na.rm = F), 
         MSL = mean(MSL, na.rm = F)) %>%
  select(Year, MHW, MSL) %>%
  mutate(TidalRange = MHW-MSL)

tideRangeLastDatum <- annapolisTideGaugeTimeSeries %>%
  filter(Year >= 1983 & Year <= 2001)

longTermTideRange <- mean(tideRangeLastDatum$TidalRange)

fit.lm2 <- lm(TidalRange~sin(2*pi*(Year-1983)/18.61) + cos(2*pi*(Year-1983)/18.61), 
              data = tideRangeLastDatum)

annapolisAmp <- summary(fit.lm2)$coefficients[1]

TargetYear <- 1920:2019

TidalRangeAtTime <- function(MHW, MSL, lunarNodalAmp, Year) {
  TidalRange <- (MHW-MSL) + (lunarNodalAmp * (sin(2*pi*(Year-1983)/18.61)))
  return(TidalRange)
}

predictionTR <- data.frame(Year = TargetYear, TidalRange = predict_TidalRange)

ggplot(data=predictionTR, aes(x=Year, y=TidalRange)) +
  geom_line(color="red")
  
  
inEqCoef <- 1
outOfEqCoef <- 1.5
rootShoot <- 2
Bmax <- 2500
# dGrowingMin <-
# dGrowingMax <- 
# tidalRange 
# SSC
# CaptureCoef
# 
