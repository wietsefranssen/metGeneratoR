#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(metGeneratoR))

if (!is.null(args)) {
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
} else {
  varType <- "tair"
  fileType <- "xy"
  nhourly <- 1
  inFile <-"~/tn_19900101.nc"
  varname <- "tn"
  if (varType == "tair") {
    inFile2 <- "~/tx_19900101.nc"
    varname2 <- "tx"
    outFile <-"~/tair_out.nc"
  } else {
    outFile <- "~/tair_out.nc"
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
}

timezone <- 0

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
# 
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
# 
# if (varType == "tair") {
#   # set_min_max_hour_c
#   hour = hour(dates)
#   minute = minute(dates)
#   yday = yday(dates)
#   dataOut <- array(NA, dim=c(nx, ny, (24 / nhourly)))
#   potrad <- array(NA, dim=c(nx, ny, (24 / nhourly)))
#   iy<-1
#   for (iy in 1:ny) {
#     ix<-1
#     lat <- lat_2d[ix,iy]
#     for (ix in 1:nx) {
#       lon <- lon_2d[ix,iy]
#       # potential_radiation_day
#       if (!is.na(indata[ix,iy])) {
#         potradtmp <- potential_radiation_day(hour, minute, yday, lon, lat, timezone)
#         # potrad[ix,iy,] <- potradtmp
#         TminHour <- (which.min(potradtmp)-1) * nhourly
#         TmaxHour <- (which.max(potradtmp)-1) * nhourly
#         dataOut[ix,iy,] <- HourlyT_cr((24 / nhourly), TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy])
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
  iy<-1
  for (iy in 1:ny) {
    ix<-1
    lat <- lat_2d[ix,iy]
    for (ix in 1:nx) {
      lon <- lon_2d[ix,iy]
      # potential_radiation_day
      if (!is.na(indata[ix,iy])) {
        # potradtmp <- potential_radiation_day(hour, minute, yday, lon, lat, timezone)
        # potrad[ix,iy,] <- potradtmp
        # TminHour <- (which.min(potradtmp)-1) * nhourly
        # TmaxHour <- (which.max(potradtmp)-1) * nhourly
        dataOut[ix,iy,] <- potential_radiation_day(hour, minute, yday, lon, lat, timezone)
      }
    }
  }
}

# iix <- 200
# iiy <- 400
# 
# plot(dataOut[iix,iiy,])
# plot(potradxx[iix,iiy,])
# plot(dataOuttmp[iix,iiy,])
# TminHourr <- (which.min(potradxx[iix,iiy,])-1)
# TmaxHourr <- (which.max(potradxx[iix,iiy,])-1)
# plot(HourlyT_cr(24, TminHourr, indata[iix,iiy], TmaxHourr, indata2[iix,iiy]))
# plot(HourlyT_cr(24, TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy]))
# t_1hourly<-HourlyT_cr(24, TminHour, indata[ix,iy], TmaxHour, indata2[ix,iy])
# t_1hourly<-HourlyT_cr(24, TminHourr, indata[iix,iiy], TmaxHourr, indata2[iix,iiy])
# plot(t_1hourly,x=c(0:23))
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
ncvar_put(ncid_out, varData, vals = dataOut)
nc_close(ncid_out)
