#!/usr/bin/env Rscript

## Load library
suppressPackageStartupMessages(library(metGeneratoR))

tminTmax <- function() {
  sunset <- 99
  noon <- 99
  sunrise <- 99
  setSunset <- F
  ## If all values are the same (eg no rediation for the whole day)
  if (length(unique(potradtmp)) == 1) {
    TminHour <- 0
    TmaxHour <- 12
  } else {
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
    
    TminHour <- sunset - 1
    TmaxHour <- floor((noon + sunrise) / 2) - 1
  }
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

potradGiven <- F
# potradFile <- "~/potrad_19900101.nc"
# potradVarname <- "potrad"

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
  datesPotrad <- seq(nc_info$ts, by = paste(1, "hours"), length = (24) * length(nc_info$ts))
  if (potradGiven) {
    potrad <- indataPotrad
    rm(indataPotrad)
  } else {
    hour = hour(datesPotrad)
    minute = minute(datesPotrad)
    yday = yday(datesPotrad)
    # potrad <- array(NA, dim=c(nx, ny, (24 / nhourly)))
    potrad <- array(NA, dim=c(nx, ny, 24))
    iy<-1
    for (iy in 1:ny) {
      ix<-1
      lat <- lat_2d[ix,iy]
      for (ix in 1:nx) {
        lon <- lon_2d[ix,iy]
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

if (varType == "tair") {
  hour = hour(dates)
  minute = minute(dates)
  yday = yday(dates)
  dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  tair_1h_arr <- array(NA, dim=c(nx, ny, (24)))
  TminHour_arr <- array(NA, dim=c(nx, ny))
  TmaxHour_arr <- array(NA, dim=c(nx, ny))
  iy<-1
  for (iy in 1:ny) {
    ix<-1
    lat <- lat_2d[ix,iy]
    for (ix in 1:nx) {
      lon <- lon_2d[ix,iy]
      if (!is.na(indata[ix,iy])) {
        potradtmp <- potrad[ix,iy,]
        res <- tminTmax()
        TminHour <- res$TminHour
        TmaxHour <- res$TmaxHour
        tair_1h <- HourlyT_cr((24), TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy])
        t_6h <- NULL
        for(i in 1:(24/nhourly)) {
          end <- i * nhourly
          start <- (end - nhourly) + 1
          dataOut[ix,iy,i] <- mean(tair_1h[start:end])
        }
        # TminHour_arr[ix,iy] <- TminHour
        # TmaxHour_arr[ix,iy] <- TmaxHour
        # tair_1h_arr[ix,iy,] <- tair_1h
      }
    }
  }
}

# iix <- 200
# iiy <- 400
# 
# lon<-lon_2d[iix,iiy]
# lat<-lat_2d[iix,iiy]
# 
# TminHour_arr[iix,iiy] 
# TmaxHour_arr[iix,iiy] 
# 
# par(mfrow=c(2,1))
# plot(potrad[iix, iiy, ], main = "potrad",x=c(0:23), xlab = "", ylab = "")
# 
# plot(tair_1h_arr[iix, iiy, ], x=c(0:23), main = "tair", xlab = "", ylab = "")
# # points(y=tair_1h_arr[iix,iiy,c(0,6,12,18)+1], x = c(0,6,12,18), pch = 20, col = "red")
# # points(y=tair_1h_arr[iix,iiy,c(5,11,17, 23)+1], x = c(5,11,17, 23), pch = 20, col = "blue")
# points(y=dataOut[iix, iiy, ], x = c(0,6,12,18), pch = 22, col = "red")
# # points(y=dataOut[iix, iiy,], x = c(5,11,17, 23), pch = 22, col = "blue")

## Add data to file
ncid_out <- nc_open(outFile, write = T )
ncvar_put(ncid_out, outVar, vals = dataOut)
nc_close(ncid_out)
