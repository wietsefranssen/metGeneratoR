rm(list = ls())
source("/home/wietse/Documents/RProjects/metGeneratoR/melo_functions.R")

library(ncdf4)
library(ncdf4.helpers)
library(units)
library(lubridate)

 # inFile <- "~/sw_1day.nc"
# varname <- "SWdown"
# fileType <- "regulairLatlon"
## cdo setgridtype,curvilinear pr_19900101_efas.nc pr_19900101_efas_curv.nc
inFile <- "~/pr_19900101_efas_curv.nc"
varname <- "pr"
fileType <- "curvilinear_2d"
# inFile <- "~/pr_19900101_efas.nc"
# varname <- "pr"
# fileType <- "xy"
outFile <- "~/Dest.nc"
timezone <- 1
nhourly <- 6

doRad <- F
doPr <- T
file.remove(outFile)


ncid_in<-nc_open(inFile, write = F)
if (fileType == "regulairLatlon") {
  inlons <- ncvar_get(ncid_in,"lon")
  nx <- length(inlons)
  #inlonsatt <- ncatt_get(ncid_in,"lon")
  inlats <- ncvar_get(ncid_in,"lat")
  ny <- length(inlats)
  #inlatsatt <- ncatt_get(ncid_in,"lat")
  ## Make d2 lon and lat arrays
  ## First based on 1d lon lat files:
  lon_1d <- inlons
  lat_1d <- inlats
  lon_2d <- array(data = NA, dim = c(nx,ny))
  lat_2d <- array(data = NA, dim = c(nx,ny))
  for (ix in 1:nx) lon_2d[ix,] <- lon_1d[ix]  
  for (iy in 1:ny) lat_2d[,iy] <- lat_1d[iy]  
} else if (fileType == "curvilinear_2d") {
  lon_2d <- inlons <- ncvar_get(ncid_in,"lon")
  nx <- ncid_in$dim$x$len
  #inlonsatt <- ncatt_get(ncid_in,"lon")
  lat_2d <- inlats <- ncvar_get(ncid_in,"lat")
  ny <- ncid_in$dim$y$len
  #inlatsatt <- ncatt_get(ncid_in,"lat")
} else {
  stop("fileType not regonised")
}

intimes <- ncvar_get(ncid_in,"time")
#intimesatt <- ncatt_get(ncid_in,"time")
ts <- nc.get.time.series(ncid_in)
indata <- ncvar_get(ncid_in, varname)
#indataatt <- ncatt_get(ncid_in,varname)
nc_close(ncid_in)


## Define target dates
dates <- seq(ts, by = paste(nhourly, "hours"), length = (24 / nhourly) * length(ts))

## Do calculation
## Radiation
if (doRad) {
  dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  iy<-1
  for (iy in 1:ny) {
    ix<-1
    lat <- lat_2d[ix,iy]
    print(paste(iy,lat))
    for (ix in 1:nx) {
      lon <- lon_2d[ix,iy]
      dataOut [ix,iy,] <- disaggregate_radiation(radiation=indata[ix,iy],date=dates,lon=lon,lat=lat,timezone=timezone)
    }
  }
}

if (doPr) {
  dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
  iy<-1
  for (iy in 1:ny) {
    ix<-1
    lat <- lat_2d[ix,iy]
    print(paste(iy,lat))
    for (ix in 1:nx) {
      lon <- lon_2d[ix,iy]
      dataOut [ix,iy,] <- indata[ix,iy]
    }
  }
}

outTimes <- seq(as.POSIXct(ts[1]), by = paste(nhourly, "hours"), length = (24 / nhourly) * length(ts))
outTatt <- paste0("hours since ", outTimes[1])
outTvals <- seq(from = 0, length.out = length(outTimes), by = nhourly)

ncid_in<-nc_open(inFile, write = F)
if (fileType == "regulairLatlon") {
  dimX <- ncdim_def(name='lon', units='degrees_east', longname='longitude', vals=inlons )
  dimY <- ncdim_def(name='lat', units='degrees_north', longname='latitude', vals=inlats )
} else if (fileType == "curvilinear_2d") {
  print("add curv")
  dimX <- ncdim_def(name='x', units='', vals = 1:nx,create_dimvar=F)
  dimY <- ncdim_def(name='y', units='', vals = 1:ny,create_dimvar=F )
  dimNV4 <- ncdim_def(name='nv4', units='', vals = 1:4,create_dimvar=F )
}
dimT <- ncdim_def(name='time', units=outTatt, vals=outTvals )

varData <- ncvar_def(name=varname, units='-', dim=list(dimX,dimY,dimT), missval=NA, prec='double')

if (fileType == "curvilinear_2d") {
  print("add curv")
  varLat <- ncvar_def(name='lat', units='degrees_north', dim=list(dimX,dimY), missval=NA, prec='double')
  varLon <- ncvar_def(name='lon', units='degrees_east', dim=list(dimX,dimY), missval=NA, prec='double')
  varLatBnds <- ncvar_def(name='lat_bnds', units='degrees_north', dim=list(dimX,dimY,dimNV4), missval=NA, prec='double')
  varLonBnds <- ncvar_def(name='lon_bnds', units='degrees_east', dim=list(dimX,dimY,dimNV4), missval=NA, prec='double')
  ncid_out <- nc_create( outFile, list(varLat,varLon,varLatBnds,varLonBnds,varData) )
} else {
  ncid_out <- nc_create( outFile, list(varData) )
}

nc_close(ncid_out)

## Add data to file
ncid_out <- nc_open( outFile, write = T )
if (fileType == "curvilinear_2d") {
  print("add curv")
  ncvar_put(ncid_out, varLon, vals = inlons)
  ncvar_put(ncid_out, varLat, vals = inlats)
  ncvar_put(ncid_out, varData, vals = dataOut)
} else {
  ncvar_put(ncid_out, varData, vals = dataOut)
}
nc_close(ncid_out)
