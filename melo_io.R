rm(list = ls())

library(ncdf4)
library(ncdf4.helpers)
library(units)
library(lubridate)

varname <- "SWdown"
destFile <- "~/Dest.nc"
nhourly <- 3


file.remove(destFile)
# ncid_in<-nc_open("~/sw_1day.nc", write = F)
ncid_in<-nc_open("~/sw_2days.nc", write = F)
inlons <- ncvar_get(ncid_in,"lon")
#inlonsatt <- ncatt_get(ncid_in,"lon")
inlats <- ncvar_get(ncid_in,"lat")
#inlatsatt <- ncatt_get(ncid_in,"lat")
intimes <- ncvar_get(ncid_in,"time")
#intimesatt <- ncatt_get(ncid_in,"time")
ts <- nc.get.time.series(ncid_in)
#indataatt <- ncatt_get(ncid_in,varname)

nc_close(ncid_in)

outTimes <- seq(as.POSIXct(ts[1]), by = paste(nhourly, "hours"), length = (24 / nhourly) * length(ts))
outTatt <- paste0("hours since ", outTimes[1])
outTvals <- seq(from = 0, length.out = length(outTimes), by = 3)

dimX <- ncdim_def(name='lon', units='degrees_east', longname='longitude', vals=inlons )
dimY <- ncdim_def(name='lat', units='degrees_north', longname='latitude', vals=inlats )
dimT <- ncdim_def(name='time', units=outTatt, vals=outTvals )

varData <- ncvar_def(name=varname, units='-', dim=list(dimX,dimY,dimT), missval=NA, prec='double')
# varLon <- ncvar_def(name='lon', units='degrees_east', dim=list(dimCross), missval=NA, longname='latitude', prec='double')

ncid_out <- nc_create( destFile, list(varData) )

nc_close(ncid_out)
