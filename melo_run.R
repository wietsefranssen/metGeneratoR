#!/usr/bin/env Rscript

## Load library
suppressPackageStartupMessages(library(metGeneratoR))

tminTmax <- function() {
  
  setSunset <- F
  for (i in 1:24) {
    i_prev <- i - 1
    if (i_prev == 0) i_prev <- 24
    pot_curr <- potradtmp[i]
    pot_prev <- potradtmp[i_prev]
    # print(paste(pot_curr, pot_prev))
    if (pot_curr < pot_prev) {
      sunrise <- i
    }
    if (pot_curr > pot_prev) {
      if (!setSunset) {
        setSunset <- T
        sunset <- i - 1
        if (sunset == 0) sunset <- 24
      }
      noon <- i;
    }
  }
  TminHour <- sunset
  TmaxHour <- floor((noon + sunrise) / 2)
  
  # print(paste("sunrise:", sunrise))
  # print(paste("sunset:", sunset))
  # print(paste0("noon: ", noon))
  result<-NULL
  result$TminHour <- TminHour
  result$TmaxHour <- TmaxHour
  return(result)
}

## Get settings
settings <- melo_get_settings()
inFile <- settings$inFile
outFile <- settings$outFile
fileType <- settings$fileType
varType <- settings$varType
nhourly <- settings$nhourly
inFile <- settings$inFile
varname <- settings$varname
if (varType == "tair") inFile2 <- settings$inFile2
if (varType == "tair") varname2 <- settings$varname2
outVar <- settings$outVar

potradGiven <- T
potradFile <- "~/potrad_19900101.nc"
potradVarname <- "potrad"

timezone <- 0
shiftinhours <- -30
shiftouthours <- 0

## Get NetCDF/Grid infomation
nc_info <- melo_get_nc_info()

## Create NetCDF
melo_create_nc(nc_info)
nx <- nc_info$nx
ny <- nc_info$ny
lon_2d <- nc_info$lon_2d
lat_2d <- nc_info$lat_2d

## Load data
ncid_in<-nc_open(inFile)
indata <- ncvar_get(ncid_in, varname)
nc_close(ncid_in)
if (varType == "tair") {
  ncid_in<-nc_open(inFile2)
  indata2 <- ncvar_get(ncid_in, varname2)
  nc_close(ncid_in)
}

## Load potrad data
if (potradGiven) {
  ncid_in<-nc_open(potradFile)
  indataPotrad <- ncvar_get(ncid_in, potradVarname)
  nc_close(ncid_in)
}

## Define dates
dates <- seq(nc_info$ts, by = paste(nhourly, "hours"), length = (24 / nhourly) * length(nc_info$ts))

## Do calculation
if (varType == "shortwave" || (varType == "tair") || varType == "potrad") {
  if (potradGiven) {
    potrad <- indataPotrad
  } else {
    hour = hour(dates)
    minute = minute(dates)
    yday = yday(dates)
    # dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
    potrad <- array(NA, dim=c(nx, ny, (24 / nhourly)))
    iy<-1
    for (iy in 1:ny) {
      ix<-1
      lat <- lat_2d[ix,iy]
      for (ix in 1:nx) {
        lon <- lon_2d[ix,iy]
        # potential_radiation_day
        if (!is.na(indata[ix,iy])) {
          potrad[ix,iy,] <- potential_radiation_day(hour, minute, yday, lon, lat, timezone)
        }
      }
    }
  }
  if (varType == "potrad") {
    dataOut <- potrad
  }
}


## Radiation
if (varType == "shortwave") {
  dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  hour = hour(dates)
  minute = minute(dates)
  yday = yday(dates)
  iy<-1
  for (iy in 1:ny) {
    ix<-1
    lat <- lat_2d[ix,iy]
    print(paste(iy,ny))
    for (ix in 1:nx) {
      lon <- lon_2d[ix,iy]
      if (!is.na(indata[ix,iy])) dataOut [ix,iy,] <- disaggregate_radiation(radiation=indata[ix,iy], hour, minute, yday, lon, lat, timezone)
    }
  }
}

if (varType == "precip") {
  dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  iy<-1
  for (iy in 1:ny) {
    ix<-1
    lat <- lat_2d[ix,iy]
    for (ix in 1:nx) {
      lon <- lon_2d[ix,iy]
      dataOut [ix,iy,] <- indata[ix,iy]
    }
  }
}

# if (varType == "tair") {
#   # set_min_max_hour_c
#   hour = 1:24
#   minute = rep(0,24)
#   yday = rep(yday(dates)[1],24)
#   dataOut <- array(NA, dim=c(nx, ny, (24)))
#   potradxx <- array(NA, dim=c(nx, ny, 24))
#   dataOuttmp <- array(NA, dim=c(nx, ny, 24))
#   iy<-1
#   for (iy in 1:ny) {
#     ix<-1
#     lat <- lat_2d[ix,iy]
#     for (ix in 1:nx) {
#       lon <- lon_2d[ix,iy]
#       # potential_radiation_day
#       if (!is.na(indata[ix,iy])) {
#         potradtmp <- potential_radiation_day(hour, minute, yday, lon, lat, timezone)
#         potradxx[ix,iy,] <- potradtmp
#         TminHour <- (which.min(potradtmp)-1) * nhourly
#         TmaxHour <- (which.max(potradtmp)-1) * nhourly
#         dataOuttmp[ix,iy,] <- HourlyT_cr(24, TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy])
#       }
#     }
#   }
# }

if (varType == "tair") {
  # set_min_max_hour_c
  hour = hour(dates)
  minute = minute(dates)
  yday = yday(dates)
  dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  # potrad <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  TminHour_arr <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  TmaxHour_arr <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  iy<-1
  for (iy in 1:ny) {
    ix<-1
    lat <- lat_2d[ix,iy]
    for (ix in 1:nx) {
      lon <- lon_2d[ix,iy]
      # potential_radiation_day
      if (!is.na(indata[ix,iy])) {
        potradtmp <- potrad[ix,iy,]
        # if (potradGiven) {
        #   potradtmp <- indataPotrad[ix,iy,]
        #   potrad[ix,iy,] <- potential_radiation_day(hour, minute, yday, lon, lat, timezone)
        # } else {
        #   potradtmp <- potential_radiation_day(hour, minute, yday, lon, lat, timezone)
        #   potrad[ix,iy,] <- potradtmp
        # }
        res<-tminTmax()
        TminHour<-res$TminHour
        TminHour_arr[ix,iy,] <- TminHour
        TmaxHour<-res$TmaxHour
        TmaxHour_arr[ix,iy,] <- TmaxHour
        # TminHour <- (which.min(potradtmp)-1) * nhourly
        # TmaxHour <- (which.max(potradtmp)-1) * nhourly
        dataOut[ix,iy,] <- HourlyT_cr((24 / nhourly), TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy])
      }
    }
  }
}


iix <- 200
iiy <- 400

lon<-lon_2d[iix,iiy]
lat<-lat_2d[iix,iiy]

plot(dataOut[iix,iiy,])
plot(potrad[iix,iiy,])
plot(indataPotrad[iix,iiy,])

TminHourr <- (which.min(potrad[iix,iiy,])-1)
TmaxHourr <- (which.max(potrad[iix,iiy,])-1)
plot(HourlyT_cr(24, TminHourr, indata[iix,iiy], TmaxHourr, indata2[iix,iiy]))
plot(HourlyT_cr(24, TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy]))
t_1hourly<-HourlyT_cr(24, TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy])
t_1hourly<-HourlyT_cr(24, TminHourr, indata[iix,iiy], TmaxHourr, indata2[iix,iiy])
plot(t_1hourly,x=c(0:23))
# t_6hourly <- NULL
# t_6hourly[1] <- mean(t_1hourly[1:6])
# t_6hourly[2] <- mean(t_1hourly[7:12])
# t_6hourly[3] <- mean(t_1hourly[13:18])
# t_6hourly[4] <- mean(t_1hourly[19:24])
# plot(t_1hourly,x=c(0:23))
# points(y=t_6hourly, x = c(0,6,12,18), pch = 22, col = "red")
# points(y=t_6hourly, x = c(5,11,17, 23), pch = 22, col = "blue")
# points(y=t_1hourly[c(0,6,12,18)+1], x = c(0,6,12,18), pch = 20, col = "red")
# points(y=t_1hourly[c(5,11,17, 23)+1], x = c(5,11,17, 23), pch = 20, col = "blue")

## Middelen VOORUIT!! (rode vierkantjes)

## Add data to file
ncid_out <- nc_open(outFile, write = T )
ncvar_put(ncid_out, outVar, vals = dataOut)
nc_close(ncid_out)
