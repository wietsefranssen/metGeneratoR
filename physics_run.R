#################################################################################
## Cleanup
rm(list=ls(all=TRUE))

setwd("/home/wietse/Documents/WORKDIRS/RProj/disaggregateR/")
source("./physics.R")

library(WFRTools)
library(lattice)

## LOAD Precip
inElev <- ncLoad(file = "./WFDEI-elevation.nc",
                 varName = "elevation")

## output dimensions
nx <- length(inElev$xyCoords$x)
ny <- length(inElev$xyCoords$y)

## Define output arrays
tiny_step_per_day <- cnst$SEC_PER_DAY / cnst$SW_RAD_DT
outRadFractions <- array(data = NA, dim = c(ny, 366, tiny_step_per_day))
outFlat_potrad  <- array(data = NA, dim = c(ny, 366))
outDaylength    <- array(data = NA, dim = c(ny, 366))
outTt_max0      <- array(data = NA, dim = c(nx, ny, 366)) # vanwege dat elke gridcell een andere elevation kan hebbben

library(doParallel)
registerDoParallel(cores=3)

ny <- 5
nx <- 2


## Do solar_geom
print(paste0("solar_geom: "))
output<-foreach(iy = 1:ny) %dopar% {
  lat <- inElev$xyCoords$y[iy]
  result <- csolar_geom(lat)
}
print(paste0("Store solar_geom:"))
for (iy in 1:ny) {
  ## ADD TO OUTPUT ARRAY
  outRadFractions[iy, , ] <- output[[iy]]$tiny_rad_fract
  outFlat_potrad[iy, ] <- output[[iy]]$flat_potrad
  outDaylength[iy, ] <- output[[iy]]$daylength
}
print(Sys.time())

print(paste0("Saving..."))
save(file = "./outRadFractions_2880.Rdata", outRadFractions)
save(file = "./outFlat_potrad_2880.Rdata", outFlat_potrad)
save(file = "./outDaylength_2880.Rdata", outDaylength)

## Do solar_geom_tt_max0
print(paste0("solar_geom_tt_max0())"))
for (iy in (1:ny) ) {
 print(paste0("iy: ", iy, ", time: " ,Sys.time()))
 output<-foreach(ix = 1:nx) %dopar% {
   if (!is.na(inElev$Data[ix,iy])) {
     lat <- inElev$xyCoords$y[iy]
     result <- csolar_geom_tt_max0(inElev$Data[ix, iy], lat, 0.0065)
   }
 }
 for (ix in 1:nx) {
   if (!is.na(inElev$Data[ix,iy])) {
     ## ADD TO OUTPUT ARRAY
     outTt_max0[ix, iy, ] <- output[[ix]]
   }
 }
}
save(file = "./outTt_max0.Rdata", outTt_max0)





# out<-csolar_geom(8.75)
# plot(out$tiny_rad_fract[60,])
plot(outRadFractions[1,,81])

# for (i in 1:80) {
#   print(plot(outRadFractions[1,,i]))
# }
outRadFractions_old<-outRadFractions

result <- csolar_geom(-89.75)
result$tiny_rad_fract[1, ]
plot(result$tiny_rad_fract[1, ])
