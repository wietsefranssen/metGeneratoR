source("/home/wietse/Documents/RProjects/metGeneratoR/melo_functions.R")

library(ncdf4)
library(ncdf4.helpers)
library(units)
library(lubridate)
library(ncmeta)
library(proj4)

# inFile <- "~/sw_1day.nc"
# varname <- "SWdown"
# fileType <- "regulairLatlon"
## cdo setgridtype,curvilinear pr_19900101_efas.nc pr_19900101_efas_curv.nc
# inFile <- "~/pr_19900101_efas_curv.nc"
# varname <- "pr"
# fileType <- "curvilinear_2d"
inFile <- "~/pr_19900101_efas.nc"
varname <- "pr"
fileType <- "xy"
outFile <- "~/Dest.nc"
timezone <- 1
nhourly <- 6

ncid_in<-nc_open(inFile, write = F)
if (fileType == "regulairLatlon") {
  inlons <- ncvar_get(ncid_in,"lon")
  nx <- length(inlons)
  inlats <- ncvar_get(ncid_in,"lat")
  ny <- length(inlats)
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
} else if (fileType == "xy") {
  inlons <- ncvar_get(ncid_in,"x")
  nx <- ncid_in$dim$x$len
  inlats <- ncvar_get(ncid_in,"y")
  ny <- ncid_in$dim$y$len
  indataatt <- ncatt_get(ncid_in,varname)
  ingridmapping <- nc_grid_mapping_atts(inFile)
  proj4string <- suppressWarnings(nc_gm_to_prj(ingridmapping))
  proj4varname <- ingridmapping$variable[1]
  x_2d <- array(data = NA, dim = c(nx,ny))
  y_2d <- array(data = NA, dim = c(nx,ny))
  for (ix in 1:nx) x_2d[ix,] <- inlons[ix]  
  for (iy in 1:ny) y_2d[,iy] <- inlats[iy]  
  x_1d <- as.vector(x_2d)
  y_1d <- as.vector(y_2d)
  xy <- data.frame(x=x_1d, y=y_1d)
  
  ## Transformed data
  pj <- proj4::project(xy, proj4string, inverse=TRUE)
  pj_lon <- pj$x
  dim(pj_lon) <- c(nx,ny)
  pj_lat <- pj$y
  dim(pj_lat) <- c(nx,ny)
  lat_2d <- pj_lat
  lon_2d <- pj_lon
} else {
  stop("fileType not regonised")
}
intimes <- ncvar_get(ncid_in,"time")
ts <- nc.get.time.series(ncid_in)
indata <- ncvar_get(ncid_in, varname)
nc_close(ncid_in)

## Define target dates
dates <- seq(ts, by = paste(nhourly, "hours"), length = (24 / nhourly) * length(ts))


outTimes <- seq(as.POSIXct(ts[1]), by = paste(nhourly, "hours"), length = (24 / nhourly) * length(ts))
outTatt <- paste0("hours since ", outTimes[1])
outTvals <- seq(from = 0, length.out = length(outTimes), by = nhourly)

ncid_in<-nc_open(inFile, write = F)
if (fileType == "regulairLatlon") {
  dimX <- ncdim_def(name='lon', units='degrees_east', longname='longitude', vals=inlons )
  dimY <- ncdim_def(name='lat', units='degrees_north', longname='latitude', vals=inlats )
} else if (fileType == "xy") {
  print("add xy")
  dimX <- ncdim_def(name='x', units='Meter', longname='x coordinate of projection', vals=inlons )
  dimY <- ncdim_def(name='y', units='Meter', longname='y coordinate of projection', vals=inlats )
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
} else if (fileType == "xy") {
  print("add xy2")
  varCoord <- ncvar_def(name=proj4varname, dim=list(), units="", prec='integer')
  ncid_out <- nc_create( outFile, list(varCoord, varData) )
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
} else if (fileType == "xy") {
  print("add xy3")

  iattNames<-names(indataatt)
  i <- 1
  for (iatt in indataatt) {
    ncatt_put(ncid_out, varname, iattNames[i], iatt)
    i <- i+1
  }
  ## add PROJ4 attributes
  for (i in 1:length(ingridmapping$name)) ncatt_put(ncid_out, proj4varname, ingridmapping$name[[i]], ingridmapping$value[[i]])
} else {
}
nc_close(ncid_out)
