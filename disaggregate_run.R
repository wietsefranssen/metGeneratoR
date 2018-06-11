## Cleanup
rm(list=ls(all=TRUE))

setwd("/home/wietse/Documents/WORKDIRS/RProj/disaggregateR/")

source("./disaggregate.R")
library(WFRTools)

cnst$SW_RAD_DT<-30
# cnst$SW_RAD_DT<-30*10

########## Calculate Rad_chunk
params <- NULL
params$time_step <- 60*3
ts <- params$time_step
ts_per_day <- (cnst$HOURS_PER_DAY * cnst$MIN_PER_HOUR / ts)

load(file = "./hpc3/outDaylength_2880.Rdata")
load(file = "./hpc3/outRadFractions_2880.Rdata")

ny <- dim(outRadFractions)[1]

rad_chunk <- array(0, dim = c(ny, 366, ts_per_day) )

for (iy in 1:ny) {
  rad_chunk[iy,,] <- tiny_rad_fract_chunk(outDaylength[iy, ], outRadFractions[iy,,], params)
}
rad_chunk[is.na(rad_chunk)] <- 0
plot(rad_chunk[,1 ,3 ])

######### Do radiation
inRsds <- ncLoad(file = "./hpc/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc", varName = "rsds", timesteps = c(1:8))
inElev <- ncLoad(file = "./WFDEI-elevation.nc", varName = "elevation")

nx <- length(inRsds$xyCoords$x)

outRsds <- inRsds
day<-1
for (iy in 1:ny) {
  for (ix in 1:nx) {
    if (!is.na(inElev$Data[ix,iy])) {
      outRsds$Data[ix,iy,] <-shortwave_day_fast(inRsds$Data[ix, iy, 1 ], rad_chunk[iy, day, ], outDaylength[iy, day])
    } else
      outRsds$Data[ix,iy,] <- NA
  }
}

ncWrite("out.nc", outRsds)

# outRsds$Data[outRsds$Data>1e18] <- NA
# image(outRsds$Data[,,5])
# image(outRsds$Data[,,5]-outRsds$Data[,,4])
#
# plot(rad_chunk[1,3,150:160])
# rad_chunk[1,1,1]
# 1:36
# 37:72
# 73:108
#
# 1:360
# 361:720
# 721:1080
#
#
# tot<-rep(NA, 360)
# for (i in 1:360) {
#   tot[i]<-sum(outRadFractions[1,1:36,i])
# }
# plot(tot)
# tot[1]
#
# source("./physics.R")
# result <- csolar_geom(-89.75)
# sum(result$tiny_rad_fract[1,1:36])
#
# result <- csolar_geom(-60.25)
# result$tiny_rad_fract[1,1:36]
#
# outRadFractions[1,1:36,1]
# csolar_geom(-89.75)$tiny_rad_fract[1,1:36]
#
# sum(outRadFractions[1,73:108,152])
# sum(csolar_geom(-14.75)$tiny_rad_fract[1,73:108])
# sum(csolar_geom(-14.25)$tiny_rad_fract[1,73:108])
# sum(csolar_geom(-13.75)$tiny_rad_fract[1,73:108])
#
# sum(csolar_geom(-14.75)$tiny_rad_fract[1,73])
# sum(csolar_geom(-14.25)$tiny_rad_fract[1,73])
# sum(csolar_geom(-13.75)$tiny_rad_fract[1,73])
#
# hh<-csolar_geom(-14.75)
# csolar_geom_old(elev = 10, lat = -61.25,lr = 0)$tiny_rad_fract[1,1:36]
#
# for (i in 1:720) {
#   csolar_geom(-61.25)$tiny_rad_fract[1,1:36]
# }

