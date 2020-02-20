rm(list = ls())
source("/home/wietse/Documents/RProjects/metGeneratoR/melo_functions_rad.R")
source("/home/wietse/Documents/RProjects/metGeneratoR/melo_functions_temp.R")

library(ncdf4)
library(ncdf4.helpers)
library(units)
library(lubridate)
library(ncmeta)
library(proj4)

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

outFile <- "~/Dest.nc"
timezone <- 1
nhourly <- 3

doRad <- T
doPr <- F

file.remove(outFile)

## Create NetCDF
source("/home/wietse/Documents/RProjects/metGeneratoR/melo_create_nc.R")

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

## Add data to file
ncid_out <- nc_open(outFile, write = T )
ncvar_put(ncid_out, varData, vals = dataOut)
nc_close(ncid_out)
