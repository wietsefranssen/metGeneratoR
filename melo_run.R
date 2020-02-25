#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(metGeneratoR))

args = commandArgs(trailingOnly=TRUE)

# Test arguments
if (length(args)<6) {
  stop("Too less arguments. \n\nRequired:\n  --> varType, projection, nhourly, inFile, varname, outfile", call.=FALSE)
} else if (args[1]=="tair") {
  if (length(args)<8) {
    stop("Too less arguments. \n\nRequired:\n  --> varType, projection, nhourly, inFileTmin, varname, inFileTmax, varname, outfile\n\n ex: ./melo_run.R tair xy 3 ~/tn_19900101.nc tn ~/tx_19900101.nc tx ~/out.nc", call.=FALSE)
  }
} else if (length(args)==2) {
  # default output file
  print(args[2])
}

varType <- args[1]
fileType <- args[2]
nhourly <- as.numeric(args[3])
inFile <- args[4]
varname <- args[5]
if (varType == "tair") {
inFile2 <- args[6]
varname2 <- args[7]
outFile <- args[8]
} else {
  outFile <- args[6]
}

## Check
if (varType != "precip" && varType != "tair" && varType != "shortwave") stop("varType not supported!")

## Print setting
cat(paste("\n######################################################################################"))
cat(paste("\n## varType:   ", varType, "    \tfileType:   ", fileType, "   nhourly:    ", nhourly))
cat(paste("\n## inFile:   ", inFile,"    varname:   ", varname))
if (varType == "tair") {
  cat(paste("\n## inFile2:   ", inFile2,"    varname2:   ", varname2))
}
cat(paste("\n## outFile:   ", outFile))
cat(paste("\n######################################################################################\n"))

timezone <- 1

## Remove the output file first (if there is one...)
if (file.exists(outFile)) {
  invisible(file.remove(outFile))
}

## Create NetCDF
source(system.file("extdata", "melo_create_nc.R",      package = "metGeneratoR"))

## Load data
ncid_in<-nc_open(inFile)
indata <- ncvar_get(ncid_in, varname)
nc_close(ncid_in)
if (varType == "tair") {
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

if (varType == "tair") {
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
