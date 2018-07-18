rm(list=ls(all=TRUE))
library(metGeneratoR)

## Settings
nx <- 720
nrec <- 24
gmt_float <- 0
sLon <- -179.75
res <- 0.5
# sLon <- -10.25


nTinyStepsPerDay <- 720



## Current location
lat <- 0.25
lon <- 0.25
# lon <- -179.75

## get ilon
ilon <- ( (lon - sLon) / res ) + 1

## Define radractions for at lat 0.25 (equator) for 1jan
radfrac <- solar_geom_cr(lat, 1, nTinyStepsPerDay)

## plot radfrac over the day
plot(radfrac)

## Check and correct gmt_float 
if(gmt_float < -12 || gmt_float > 12) stop("cannot be lower than -12 and higher than 12")
if (gmt_float < 0) gmt_float <- gmt_float + nTinyStepsPerDay

## Define the index offset based on the gmt offset
gmt_float_tmp <- gmt_float + (24/nrec)/2
iGmtOffset <- gmt_float_tmp * (nTinyStepsPerDay/24)


## Define the index offset longitude location
iLonOffset <- (( (lon / res) + res ) )/ (720/nTinyStepsPerDay)
if (iLonOffset < 0) iLonOffset <- iLonOffset + nTinyStepsPerDay
iLonOffset

radfrac_new <- array(NA, dim = nTinyStepsPerDay)
radfrac_new_rec <- array(0, dim = nrec)
for (i in 1:nTinyStepsPerDay) {
  iNew <- i + iGmtOffset + iLonOffset
  if (iNew > nTinyStepsPerDay) iNew <- (iNew - nTinyStepsPerDay)
  radfrac_new[iNew] <- radfrac[i]
}

for (i in 1:nTinyStepsPerDay) {
  ## rec...
  rec <- ceiling(i/ (nTinyStepsPerDay / nrec))
  radfrac_new_rec[rec] <- radfrac_new_rec[rec] + radfrac_new[i]
}

plot(radfrac)
plot(radfrac_new)
plot(radfrac_new_rec)
max(radfrac_new_rec)
