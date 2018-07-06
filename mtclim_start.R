rm(list=ls(all=TRUE))

library(metGeneratoR)
profile<-NULL

mgsetLonlatbox(c(92.25, 110.25, 7.25, 36.25))
# mgsetLonlatbox(c(92.25, 92.25, -8.25, -8.25))
# mgsetLonlatbox(c(92.25, 92.75, 34.25, 36.75))
# mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))
# mgsetPeriod(startdate = "1950-01-03", enddate = "1950-01-05")
# mgsetPeriod(startdate = "1964-12-31", enddate = "1965-01-3")
mgsetPeriod(startdate = "1965-01-01", enddate = "1965-01-2")
mgsetInDt(24) # Set N hours per timestep
mgsetOutDt(6) # Set N hours per timestep

metGen$constants<-setConstants()
constants <- metGen$constants

have_dewpt <- F
theta_l = -30

param_set_TYPE_SHORTWAVE_SUPPLIED <- 1
options <- NULL
options$SW_PREC_THRESH <- 0
options$MTCLIM_SWE_CORR <- 0
options$LW_CLOUD <- 1
options$VP_ITER <- 1

mgsetInVars(list(
  pr         = list(ncname = "pr",     filename = "../example_data4mtclim/Global/pr_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  tasmin     = list(ncname = "tasmin", filename = "../example_data4mtclim/Global/tasmin_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  tasmax     = list(ncname = "tasmax", filename = "../example_data4mtclim/Global/tasmax_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  shortwave  = list(ncname = "rsds",   filename = "../example_data4mtclim/Global/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  longwave   = list(ncname = "rlds",   filename = "../example_data4mtclim/Global/rlds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  wind       = list(ncname = "sfcWind",   filename = "../example_data4mtclim/Global/wind_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc")
))

## Define elevation file
mgsetElevation(ncname = "elevation", filename = metGen$internal$ncFileNameElevation)

## Define output variables
# mgsetOutVars(c("pr", "tas"))
mgsetOutVars(c( "shortwave", "tas"))


## Load elevation
elev <- ncLoad(file = metGen$internal$ncFileNameElevation,
               varName = metGen$settings$elevation$ncname,
               lonlatbox = metGen$settings$lonlatbox)

## Load mask (or base it on elevation)
mask<-elev
mask$Data[!is.na(mask$Data)]<- 1

nx <- length(mask$xyCoords$x)
ny <- length(mask$xyCoords$y)

## Calculate solar GEOMs as preprocessing step
lapse_rate <- 0.0065
solar_geom_file <- "./solar_geom_new.Rdata"
if(file.exists(solar_geom_file)) {
  printf("Loading solar_geom_file: %s\n", solar_geom_file)
  load(file = solar_geom_file)
} else {
  print("No solar_geom_file file!...") 
}
## Check the file:
if (!is.null(solar_geom) && dim(solar_geom$lats) == dim(mask$xyCoords$y) && dim(solar_geom$lons) == dim(mask$xyCoords$x) && dim(solar_geom$elev) == dim(mask$Data)) {
  print("File is oke!")
} else {
  print("Error in file! Regenerating solar_geom file...") 
  # solar_geom<-solar_geom(lat = lat)
  solar_geom <- NULL
  solar_geom$tiny_rad_fract <- array(NA, dim = c(366, metGen$derived$nSwRadDay, length(mask$xyCoords$y)))
  solar_geom$daylength <- array(NA, dim = c(366, length(mask$xyCoords$y)))
  solar_geom$flat_potrad <- array(NA, dim = c(366, length(mask$xyCoords$y)))
  solar_geom$tt_max <- array(NA, dim = c(366, length(mask$xyCoords$x), length(mask$xyCoords$y)))
  solar_geom$lons <- mask$xyCoords$x
  solar_geom$lats <- mask$xyCoords$y
  solar_geom$elev <- elev$Data
  
  for (ilat in 1:length(mask$xyCoords$y)) {
    print(ilat)
    solar_geom_temp <- solar_geom(lat = mask$xyCoords$y[ilat])
    solar_geom$tiny_rad_fract[,,ilat] <- solar_geom_temp$tiny_rad_fract
    solar_geom$daylength[,ilat] <- solar_geom_temp$daylength
    solar_geom$flat_potrad[,ilat] <- solar_geom_temp$flat_potrad
    for (ilon in 1:length(mask$xyCoords$x)) {
      if (!is.na(mask$Data[ilon,ilat])) {
        solar_geom$tt_max[,ilon,ilat] <- solar_geom_tt_max0(lat = mask$xyCoords$y[ilat], elev = elev$Data[ilon,ilat], lr = lapse_rate)
      }
    }
  }
  rm(solar_geom_temp)
  save(solar_geom, file = solar_geom_file)
}

### TODO: check if start yday is correct!!!
hourly_rad_fract_file <- "./hourly_rad_fract_new.Rdata"
if(file.exists(hourly_rad_fract_file)) {
  printf("Loading hourly_rad_fract_file: %s\n", hourly_rad_fract_file)
  load(file = hourly_rad_fract_file)
} else {
  print("No hourly_rad_fract_file file!...") 
}
## Check the file:
if (!is.null(hourly_rad_fract) && dim(hourly_rad_fract$lat) == dim(mask$xyCoords$y) && dim(hourly_rad_fract$lon) == dim(mask$xyCoords$x)) {
  print("File is oke!")
  hourly_rad_fract_data <- hourly_rad_fract$data
} else {
  print("Error in file! Regenerating hourly_rad_fract file...") 
  hourly_rad_fract<-NULL
  hourly_rad_fract$data <- hourly_rad_fract(mask, theta_l, solar_geom$tiny_rad_fract)
  hourly_rad_fract$lons <- mask$xyCoords$x
  hourly_rad_fract$lats <- mask$xyCoords$y
  save(hourly_rad_fract, file = hourly_rad_fract_file)
  hourly_rad_fract_data <- hourly_rad_fract$data
}

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
  # iday<-2
  metGen$current$timestep <- 1
  printf("Day: %d\n", iday)
  
  ## Init progressbar
  pb <- txtProgressBar(min = 0, max = length(mask$xyCoords$x), initial = 0, char = ">",
                       width = 80, title, label, style = 1, file = "")
  
  ## LOAD WHOLE DOMAIN FROM NETCDF
  profile$start.time.read <- Sys.time()
  inData <- readAllForcing(mask, iday)
  profile$end.time.read <- Sys.time()
  
  yday            <- metGen$derived$inYDays[iday]
  
  # for (ilat in 1:length(mask$xyCoords$y)) {
  for (ilon in 1:length(mask$xyCoords$x)) {
    for (ilat in 1:length(mask$xyCoords$y)) {
      ## Select mt data for current day in year
      # mt$ttmax0       <- solar_geom$tt_max[yday, ilon, ilat]
      flat_potrad  <- solar_geom$flat_potrad[yday, ilat]
      slope_potrad <- solar_geom$flat_potrad[yday, ilat]
      daylength    <- solar_geom$daylength[yday, ilat]
      
      elevation <- mask$Data[ilon, ilat]
      if (!is.na(elevation)) {
        # printf("Day:  %d, elevation %d, lon %5.2f, lat %5.2f\n", iday, elevation, lon, lat)
        
        lon <- mask$xyCoords$x[ilon]
        lat <- mask$xyCoords$y[ilat]
        prec <- inData$pr[ilon, ilat,1]
        tmin <- inData$tasmin[ilon, ilat,1]
        tmax <- inData$tasmax[ilon, ilat,1]
        rsds <- inData$shortwave[ilon, ilat,1]
        
        ## Calculate offset longitude
        hour_offset = (theta_l-lon)*24/360;
        if (hour_offset < 0) {
          hour_offset_int = floor(hour_offset-0.5);
        } else {
          hour_offset_int = floor(hour_offset+0.5);
        }
        hour_offset <- hour_offset - hour_offset_int #// hour_offset is now the distance from the center of local time zone
        
        # /*************************************************
        #   Shortwave, part 1.
        # *************************************************/
        hourly_rad <- array(NA, dim = c(24))
        if (param_set_TYPE_SHORTWAVE_SUPPLIED) {
          if(metGen$derived$inDt == 24) {
            hourly_rad <- rep(rsds, 24)
          } else {
            for (hour in 1:24) {
              hourly_rad[hour] = rsds[hour];
            }
          }
        }
        ###################### START MTCLIM WRAPPER
        mt <- mtclim_init(have_dewpt, param_set_TYPE_SHORTWAVE_SUPPLIED, elevation, 0,0,0,
                          lat, prec, tmax, tmin, vp, metGen$derived$inYDays[iday], hourly_rad,
                          p, mtclim_data)
        
        mt<-calc_tair(mt)
        
        mt<-calc_prcp(mt)
        
        mt<-snowpack(mt)
        
        ## Select mt data for current day in year
        mt$ttmax0       <- solar_geom$tt_max[yday, ilon, ilat]
        mt$flat_potrad  <- flat_potrad
        mt$slope_potrad <- slope_potrad
        mt$daylength    <- daylength
        
        mtclim_data<-calc_rest_1t(mt, options)
        ctrl<-mt$ctrl
        
        if (ctrl$insw) {
          tmp_rad <- mtclim_data$s_srad * 24.;
        } else {
          tmp_rad = mtclim_data$s_srad * mtclim_data$s_dayl / 3600.;
        }
        hourly_rad <- array(0, dim = (24))
        
        ## Previous mtclim_to_vic function now splitted up in hourly_rad_fract function and the lines below:
        for (j in 1:24) {
          hourly_rad[j] <-  hourly_rad_fract_data[yday,j, ilon, ilat] * tmp_rad;
        }
        
        # print(hourly_rad)
        ###################### END MTCLIM WRAPPER
        
        # /***********************************************************
        #   Shortwave, part 2.
        # Transfer the hourly shortwave from MTCLIM to atmos array.
        # This hourly shortwave is one of the following:
        #   a) exactly equal to the supplied shortwave, if supplied shortwave was hourly
        # b) equal to the supplied shortwave when aggregated up to the DT of the supplied shortwave (with hourly variability estimated by MTCLIM)
        # c) completely estimated by MTCLIM, if no shortwave was supplied as a forcing
        # ***********************************************************/
        #
        # // Ignore MTCLIM estimates if sub-daily SW was supplied
        if (param_set_TYPE_SHORTWAVE_SUPPLIED && (metGen$derived$inDt < 24)) {
          for (day in 1:metGen$derived$nday) {
            for (hour in 1:24) {
              hourly_rad[(day-1)*24+hour] = local_forcing_data[SHORTWAVE][(day-1)*24+hour];
            }
          }
        }
        
        ## Copy data to previous timestep for first timestep
        if (iday == 1) {
          hourly_rad_prev <- hourly_rad 
        }
        
        ## Fill the arrays
        # // Transfer hourly_rad to atmos structure
        for(rec in 1:metGen$derived$nOutStepDay) {
          sum = 0;
          hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
          if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
          outData$shortwave[ilon, ilat, rec] <- 0;
          for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
            if (idx_tmp <= 24) {
              idx<-idx_tmp
              outData$shortwave[ilon, ilat, rec] <- outData$shortwave[ilon, ilat, rec] + hourly_rad_prev[idx];
              # printf("rec %d, hour_local_out %d, prev_day index: %d\n", rec, hour_local_out, idx)
            } else {
              idx<-idx_tmp-24
              outData$shortwave[ilon, ilat, rec] <- outData$shortwave[ilon, ilat, rec] + hourly_rad[idx];
              # printf("rec %d, hour_local_out %d, curr_day index: %d\n", rec, hour_local_out, idx)
            }
          }
          outData$shortwave[ilon, ilat, rec] <- outData$shortwave[ilon, ilat, rec] / metGen$derived$outDt
          sum <- sum + outData$shortwave[ilon, ilat, rec]
        }
        
        ## Temperature!!
        tminmaxhour<- set_t_minmax_hour(hourly_rad)
        
        hourly_tair<-HourlyT(tminmaxhour[2],tmax,tminmaxhour[1],tmin)
        
        if (iday == 1) {
          hourly_tair_prev <- hourly_tair
        }
        
        # // Transfer hourly_tair to atmos structure
        for(rec in 1:metGen$derived$nOutStepDay) {
          sum = 0;
          hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
          if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
          outData$tas[ilon, ilat, rec] <- 0;
          for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
            if (idx_tmp <= 24) {
              idx<-idx_tmp
              outData$tas[ilon, ilat, rec] <- outData$tas[ilon, ilat, rec] + hourly_tair_prev[idx];
              # printf("rec %d, hour_local_out %d, prev_day index: %d\n", rec, hour_local_out, idx)
            } else {
              idx<-idx_tmp-24
              outData$tas[ilon, ilat, rec] <- outData$tas[ilon, ilat, rec] + hourly_tair[idx];
              # printf("rec %d, hour_local_out %d, curr_day index: %d\n", rec, hour_local_out, idx)
            }
          }
          outData$tas[ilon, ilat, rec] <- outData$tas[ilon, ilat, rec] / metGen$derived$outDt
          sum <- sum + outData$tas[ilon, ilat, rec]
        }
        
        
        ## Move data to previous timestep
        hourly_rad_prev <- hourly_rad 
        hourly_tair_prev <- hourly_tair 
        
      }
    }
    ## refresh progressbar
    setTxtProgressBar(pb, ilon)
  }
  
  ## Close ProgressBar
  close(pb)
  
  ## ADD OUTPUT TO NETCDF
  profile$start.time.write <- Sys.time()
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
  profile$end.time.write <- Sys.time()
  
}

profile$end.time.total <- Sys.time()
cat(sprintf("  Total: %.1f seconds\n", as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")))

ncellsTotal <- dim(mask$Data)[1] * dim(mask$Data)[2]
ncellsTotal <- sum(mask$Data, na.rm = TRUE)
totTime<-as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")
cat(sprintf("  Total: %.1f seconds for %d x %d = %d cells. So: 100 years (67420 cels) will take: %.1f days\n", totTime, metGen$derived$nday, ncellsTotal, 
            (metGen$derived$nday * ncellsTotal), ((totTime * 67420 * 365.25 * 100) / (metGen$derived$nday * ncellsTotal) / 86400
            )))