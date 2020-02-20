rm(list = ls())
#source("melo_functions_rad.R")
#source("/home/wietse/Documents/RProjects/metGeneratoR/melo_functions_rad.R")
# source("/home/wietse/Documents/RProjects/metGeneratoR/melo_functions_temp.R")

library(metGeneratoR)
library(ncdf4)
library(ncdf4.helpers)
library(units)
library(lubridate)
library(ncmeta)
library(proj4)

source(system.file("Rscripts", "melo_functions_rad.R",      package = "metGeneratoR"))

inFile <- "~/sw_1day.nc"
varname <- "SWdown"
fileType <- "latlon"

## cdo setgridtype,curvilinear pr_19900101_efas.nc pr_19900101_efas_curv.nc
# inFile <- "~/pr_19900101_efas_curv.nc"
# varname <- "pr"
# fileType <- "curvilinear_2d"

# inFile <- "~/pr_19900101_efas.nc"
# varname <- "pr"
# fileType <- "xy"

inFile <- "~/tmin_daily.nc"
varname <- "Tair"
inFile2 <- "~/tmax_daily.nc"
varname2 <- "Tair"
fileType <- "latlon"

inFile <- "~/tn_19920628.nc"
varname <- "tn"
inFile2 <- "~/tx_19920628.nc"
varname2 <- "tx"
fileType <- "xy"

outFile <- "~/Dest.nc"
timezone <- 1
nhourly <- 3

doRad <- F
doPr <- F
doTemp <- T

TminHour <- 7
TmaxHour <- 14


file.remove(outFile)

## Create NetCDF
# source("/home/wietse/Documents/RProjects/metGeneratoR/R/melo_create_nc.R")
source(system.file("Rscripts", "melo_create_nc.R",      package = "metGeneratoR"))

## Load data
ncid_in<-nc_open(inFile)
indata <- ncvar_get(ncid_in, varname)
nc_close(ncid_in)
if (doTemp) {
  ncid_in<-nc_open(inFile2)
  indata2 <- ncvar_get(ncid_in, varname2)
  nc_close(ncid_in)
}


## Define target dates
dates <- seq(ts, by = paste(nhourly, "hours"), length = (24 / nhourly) * length(ts))

## Do calculation
## Radiation
if (doRad) {
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

if (doPr) {
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

if (doTemp) {
  # set_min_max_hour_c
  hour = hour(dates)
  minute = minute(dates)
  yday = yday(dates)
  dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  iy<-1
  for (iy in 1:ny) {
    ix<-1
    lat <- lat_2d[ix,iy]
    for (ix in 1:nx) {
      lon <- lon_2d[ix,iy]
      # potential_radiation_day
      if (!is.na(indata[ix,iy])) {
        potradtmp <- potential_radiation_day(hour, minute, yday, lon, lat, timezone)
        TminHour <- (which.min(potradtmp)-1) * nhourly
        TmaxHour <- (which.max(potradtmp)-1) * nhourly
        dataOut [ix,iy,] <- HourlyT_cr((24 / nhourly), TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy])
      }
    }
  }
}

## Add data to file
ncid_out <- nc_open(outFile, write = T )
ncvar_put(ncid_out, varData, vals = dataOut)
nc_close(ncid_out)
