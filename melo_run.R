#!/usr/bin/env Rscript
#rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<6) {
  stop("Too less arguments. \n\nRequired:\n  --> vartype, projection, nhourly, inFile, varname, outfile", call.=FALSE)
} else if (args[1]=="temp") {
  if (length(args)<8) {
    stop("Too less arguments. \n\nRequired:\n  --> vartype, projection, nhourly, inFileTmin, varname, inFileTmax, varname, outfile\n\n ex: ./melo_run.R temp xy 3 ~/tn_19900101.nc tn ~/tx_19900101.nc tx ~/out.nc", call.=FALSE)
  }
} else if (length(args)==2) {
  # default output file
  print(args[2])
}

library(metGeneratoR)
library(ncdf4)
library(ncdf4.helpers)
library(units)
library(lubridate)
library(ncmeta)
library(proj4)
#source("/home/wietse/Documents/RProjects/metGeneratoR/melo_functions_rad.R")
source(system.file("extdata", "melo_functions_rad.R",      package = "metGeneratoR"))

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

vartype <- args[1]
fileType <- args[2]
nhourly <- args[3]
inFile <- args[4]
varname <- args[5]
inFile2 <- args[6]
varname2 <- args[7]
outFile <- args[8]

timezone <- 1
nhourly <- 3

varType <- "shortwave"
varType <- "precip"
varType <- "temp"

TminHour <- 7
TmaxHour <- 14


file.remove(outFile)

## Create NetCDF
# source("/home/wietse/Documents/RProjects/metGeneratoR/melo_create_nc.R")
source(system.file("extdata", "melo_create_nc.R",      package = "metGeneratoR"))

## Load data
ncid_in<-nc_open(inFile)
indata <- ncvar_get(ncid_in, varname)
nc_close(ncid_in)
if (varType == "temp") {
  ncid_in<-nc_open(inFile2)
  indata2 <- ncvar_get(ncid_in, varname2)
  nc_close(ncid_in)
}


## Define target dates
dates <- seq(ts, by = paste(nhourly, "hours"), length = (24 / nhourly) * length(ts))

## Do calculation
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

if (varType == "temp") {
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
