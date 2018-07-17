## ISSUES/TODO
## if eg prev_pr is used than the map is moved one gridcell upwards 
## vapor pressure is too high???
## shortwave can only be done globally: fix this!
## check if start yday is correct
## unitconversion


rm(list=ls(all=TRUE))

library(metGeneratoR)

profile<-NULL

# mgsetLonlatbox(c(92.25, 110.25, 7.25, 36.25))
# mgsetLonlatbox(c(92.25, 92.25, -8.25, -8.25))
# mgsetLonlatbox(c(92.25, 92.75, 34.25, 36.75))
mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))
mgsetPeriod(startdate = "1950-4-01", enddate = "1950-4-2")
# mgsetPeriod(startdate = "1964-12-31", enddate = "1965-01-3")
# mgsetPeriod(startdate = "1965-01-01", enddate = "1965-06-2")
mgsetInDt(24) # Set N hours per timestep
mgsetOutDt(1) # Set N hours per timestep

metGen$constants<-setConstants()
constants <- metGen$constants

theta_l = -30

options <- NULL
options$SW_PREC_THRESH <- 0
options$MTCLIM_SWE_CORR <- 0
options$LW_CLOUD <- 1
options$VP_ITER <- 1

mgsetInVars(list(
  pr         = list(ncname = "pr",      filename = "../example_data4mtclim/Global/pr_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  tasmin     = list(ncname = "tasmin",  filename = "../example_data4mtclim/Global/tasmin_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  tasmax     = list(ncname = "tasmax",  filename = "../example_data4mtclim/Global/tasmax_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  pressure   = list(ncname = "ps",      filename = "../example_data4mtclim/Global/ps_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  relhum     = list(ncname = "hurs",    filename = "../example_data4mtclim/Global/hurs_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  shortwave  = list(ncname = "rsds",    filename = "../example_data4mtclim/Global/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  longwave   = list(ncname = "rlds",    filename = "../example_data4mtclim/Global/rlds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  wind       = list(ncname = "sfcWind", filename = "../example_data4mtclim/Global/wind_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc")
))

## Define elevation file
mgsetElevation(ncname = "elevation", filename = metGen$internal$ncFileNameElevation)

## Define output variables
# mgsetOutVars(c("pr", "tas"))
mgsetOutVars(c( "shortwave", "longwave", "tas", "pr", "pressure", "wind", "vp"))
# mgsetOutVars(c( "shortwave", "pr"))

## Load elevation
elev <- ncLoad(file = metGen$internal$ncFileNameElevation,
               varName = metGen$settings$elevation$ncname,
               lonlatbox = metGen$settings$lonlatbox)

## Load mask (or base it on elevation)
mask<-elev
mask$Data[!is.na(mask$Data)]<- 1

nx <- length(mask$xyCoords$x)
ny <- length(mask$xyCoords$y)

# load("./hpc/outRadFractions_2880.Rdata")
# load("./hpc/outDaylength_2880.Rdata")
# load("./hpc/outFlat_potrad_2880.Rdata")
# load("./hpc/outTt_max0.Rdata")

## makeOutputNetCDF
makeNetcdfOut(mask)

## DEFINE OUTPUT ARRAY
outData <- NULL
for (var in names(metGen$settings$outVars)) {
  outData[[var]] <- array(NA, dim = c(nx, ny, metGen$derived$nOutStepDay))
}

### THE MAIN LOOP
profile$start.time.total <- Sys.time()
for (iday in 1:metGen$derived$nday) {
  metGen$current$timestep <- iday
  yday            <- metGen$derived$inYDays[iday]
  printf("Day: %d, date: %s\n", iday, metGen$derived$inDates[iday])
  
  ## LOAD WHOLE DOMAIN FROM NETCDF
  profile$start.time.read <- Sys.time()
  inData <- readAllForcing(mask, yday)
  inData$shortwave[,,1] <- mask$Data * inData$shortwave[,,1]
  profile$end.time.read <- Sys.time()
  
  radfrac<-aperm(rad_map_final_cr(metGen$derived$nOutStepDay, yday, nx_parts = 720), c(3,2,1))
  ## Mask out
  for (i in 1:metGen$derived$nOutStepDay) {
    ccc<-radfrac[,,i]
    radfrac[,,i] <- ccc
  }
  #image(radfrac[,,3])
  #plot(radfrac[1,,2])
  
  # /*************************************************
  #   Precipitation
  # *************************************************/
  if(!is.null(metGen$settings$inVar$pr) && !is.null(outData$pr)) {
    for(rec in 1:metGen$derived$nOutStepDay) {
      outData$pr[, , rec] <- inData$pr[, ,1]
    }
  }

  # /*************************************************
  #   Shortwave radiation
  # *************************************************/
  if(!is.null(metGen$settings$inVar$shortwave) && !is.null(outData$shortwave)) {
    for(rec in 1:metGen$derived$nOutStepDay) {
      outData$shortwave[, , rec] <- radfrac[ , , rec] * inData$shortwave[, ,1]
    }
  }
  
  # /*************************************************
  #   Longwave radiation
  # *************************************************/
  if(!is.null(metGen$settings$inVar$longwave) && !is.null(outData$longwave)) {
    for(rec in 1:metGen$derived$nOutStepDay) {
      outData$longwave[, , rec] <- radfrac[ , , rec] * inData$longwave[, ,1]
    }
  }
  
  # /*************************************************
  #   Wind
  # *************************************************/
  if(!is.null(metGen$settings$inVar$wind) && !is.null(outData$wind)) {
    for(rec in 1:metGen$derived$nOutStepDay) {
      outData$wind[, , rec] <- inData$wind[, ,1]
    }
  }
  
  # /**************************************
  #   Atmospheric Pressure (Pa)
  # **************************************/
  if(!is.null(metGen$settings$inVar$pressure) && !is.null(outData$pressure)) {
    for(rec in 1:metGen$derived$nOutStepDay) {
      outData$pressure[, , rec] <- inData$pressure[, ,1]
    }
  }
  
  # /**************************************
  #   Temperature
  # **************************************/
  if(!is.null(metGen$settings$inVar$tasmin) && !is.null(metGen$settings$inVar$tasmin) && !is.null(outData$tas)) {
    for(rec in 1:metGen$derived$nOutStepDay) {
      # tminmaxhour<- set_t_minmax_hour(hourly_rad)
      # subdaily_tair<-HourlyT(tminmaxhour[2],tmax,tminmaxhour[1],tmin)
      # 
      # outData$pressure[, , rec] <- subdaily_tair[, ,rec]
    }
  }

  # # for (ilat in 1:length(mask$xyCoords$y)) {
  # for (ilon in 1:length(mask$xyCoords$x)) {
  #   lon      <- mask$xyCoords$x[ilon]
  #   
  #   ## Calculate offset longitude
  #   nrOffsetSteps <- 24
  #   # nrOffsetSteps <- 720
  #   
  #   hour_offset <- hour_offset_int <- ceiling(ilon * (nrOffsetSteps/720))    # hour_offset<-0
  #   
  #   for (ilat in 1:length(mask$xyCoords$y)) {
  #     ## Select mt data for current day in year
  #     # mt$ttmax0       <- solar_geom$tt_max[yday, ilon, ilat]
  #     # flat_potrad  <- solar_geom$flat_potrad[yday, ilat]
  #     # slope_potrad <- solar_geom$flat_potrad[yday, ilat]
  #     # daylength    <- solar_geom$daylength[yday, ilat]
  #     
  #     elevation <- mask$Data[ilon, ilat]
  #     if (!is.na(elevation) && !is.na(inData$pr[ilon, ilat,1])) {
  #       
  #       lat      <- mask$xyCoords$y[ilat]
  #       prec     <- inData$pr[ilon, ilat,1]
  #       pressure <- inData$pressure[ilon, ilat,1]
  #       relhum   <- inData$relhum[ilon, ilat,1] / 100 # convert to fraction
  #       wind <- inData$wind[ilon, ilat,1]
  #       tmin <- inData$tasmin[ilon, ilat,1] - 273.15
  #       tmax <- inData$tasmax[ilon, ilat,1] - 273.15
  #       shortwave <- inData$shortwave[ilon, ilat,1]
  #       longwave <- inData$longwave[ilon, ilat,1]
  #       
  #       
  #       ## Temperature!!
  #       if(!is.null(metGen$settings$inVar$tasmin) && !is.null(metGen$settings$inVar$tasmin) && !is.null(outData$tas)) {
  #         tminmaxhour<- set_t_minmax_hour(hourly_rad)
  #         hourly_tair<-HourlyT(tminmaxhour[2],tmax,tminmaxhour[1],tmin)
  #         
  #         if (iday == 1) {
  #           hourly_tair_prev <- hourly_tair
  #         }
  #         
  #         # // Transfer hourly_tair to atmos structure
  #         for(rec in 1:metGen$derived$nOutStepDay) {
  #           sum = 0;
  #           hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
  #           if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
  #           outData$tas[ilon, ilat, rec] <- 0;
  #           for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
  #             if (idx_tmp <= 24) {
  #               idx<-idx_tmp
  #               outData$tas[ilon, ilat, rec] <- outData$tas[ilon, ilat, rec] + hourly_tair_prev[idx];
  #             } else {
  #               idx<-idx_tmp-24
  #               outData$tas[ilon, ilat, rec] <- outData$tas[ilon, ilat, rec] + hourly_tair[idx];
  #             }
  #           }
  #           outData$tas[ilon, ilat, rec] <- outData$tas[ilon, ilat, rec] / metGen$derived$outDt
  #           sum <- sum + outData$tas[ilon, ilat, rec]
  #         }
  #       }
  #       
  #       
  #       # /*************************************************
  #       #   Vapor pressure
  #       # *************************************************/
  #       if(!is.null(metGen$settings$inVar$relhum) && !is.null(outData$vp)) {
  #         
  #         if (iday == 1) {
  #           relhum_prev <- relhum
  #         }
  #         
  #         if(!is.null(metGen$settings$inVar$relhum) && !is.null(metGen$settings$outVars$tas)) {
  #           for(rec in 1:metGen$derived$nOutStepDay) {
  #             hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
  #             if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
  #             outData$vp[ilon, ilat, rec] <- 0;
  #             for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
  #               if (idx_tmp <= 24) {
  #                 outData$vp[ilon, ilat, rec] <- relhum_prev * svp(outData$tas[ilon, ilat, rec]) / 100
  #               } else {
  #                 outData$vp[ilon, ilat, rec] <- relhum * svp(outData$tas[ilon, ilat, rec]) / 100
  #               }
  #             }
  #           }
  #         }
  #       }
  #       
  #       ## Move data to previous timestep
  #       # hourly_rad_prev <- hourly_rad 
  #       # hourly_tair_prev <- hourly_tair 
  #       prec_prev <- prec 
  #       # pressure_prev <- pressure
  #       # wind_prev <- wind
  #       # longwave_prev <- longwave
  #       # relhum_prev <- relhum
  #     }
  #   }
  
  ## ADD OUTPUT TO NETCDF
  for (var in names(metGen$settings$outVars)) {
    # var<-"shortwave"
    timeIndex <- metGen$derived$nOutStepDay*(iday-1)+1
    metGen$settings$outVars[[var]]$ncid <- nc_open(metGen$settings$outVars[[var]]$filename, write = TRUE)
    ncvar_put(metGen$settings$outVars[[var]]$ncid,
              var,
              outData[[var]][,,],
              start = c(1, 1, timeIndex),
              count = c(nx, ny, metGen$derived$nOutStepDay)
    )
    nc_close(metGen$settings$outVars[[var]]$ncid)
  }
}

ncellsTotal <- dim(mask$Data)[1] * dim(mask$Data)[2]
ncellsTotal <- sum(mask$Data, na.rm = TRUE)
profile$end.time.total <- Sys.time()
totTime<-as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")
cat(sprintf("  Total: %.1f seconds for %d x %d = %d cells. So: 100 years (67420 cels) will take: %.1f days\n", totTime, metGen$derived$nday, ncellsTotal, 
            (metGen$derived$nday * ncellsTotal), ((totTime * 67420 * 365.25 * 100) / (metGen$derived$nday * ncellsTotal) / 86400
            )))

