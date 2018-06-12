## Cleanup
rm(list=ls(all=TRUE))

library(disaggregateR)

params <- NULL
params$time_step <- 60*3
params$do_rad_small <- F

cnst$SW_RAD_DT <- 30
# cnst$SW_RAD_DT <- 30*10

######### Do load rad_small
if(params$do_rad_small) {
  load(file = "./hpc/outDaylength_2880.Rdata")
  load(file = "./hpc/outRadFractions_2880.Rdata")
  
  rad_small <- rad_fract_chunk_small(outDaylength, outRadFractions, params)
  save(rad_small , file = "rad_small.Rdata")
  rm(outRadFractions, outDaylength)
} else {
  load(file = "rad_small.Rdata")
  load(file = "./hpc/outDaylength_2880.Rdata")
}

######### Do radiation
inRsds <- ncLoad(file = "./hpc/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc", varName = "rsds", timesteps = c(1:8))
# inRsds <- ncLoad(file = "./hpc/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc", varName = "rsds")
inElev <- ncLoad(file = "./WFDEI-elevation.nc", varName = "elevation")

nx <- length(inRsds$xyCoords$x)
ny <- length(inRsds$xyCoords$y)

outRsds <- inRsds
day<-1
for (iy in 1:ny) {
  for (ix in 1:nx) {
    if (!is.na(inElev$Data[ix,iy])) {
      outRsds$Data[ix,iy,] <- shortwave_day_fast(inRsds$Data[ix, iy, 1 ], rad_small[iy, day, ], outDaylength[iy, day])
    } else
      outRsds$Data[ix,iy,] <- NA
  }
}

ncWrite("out.nc", outRsds)

######################
settings <- mtclim_getSettings()
mtclim_run(settings)
