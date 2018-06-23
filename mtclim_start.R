rm(list=ls(all=TRUE))

library(metGeneratoR)
profile<-NULL

# mgsetLonlatbox(c(92.25, 110.25, 7.25, 36.25))
# mgsetLonlatbox(c(92.25, 92.25, -8.25, -8.25))
mgsetLonlatbox(c(92.25, 92.75, 34.25, 36.75))
# mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))
# mgsetPeriod(startdate = "1950-01-03", enddate = "1950-01-05")
# mgsetPeriod(startdate = "1964-12-31", enddate = "1965-01-3")
mgsetPeriod(startdate = "1965-01-01", enddate = "1965-01-20")
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

prec <- c(0.79, 0.00, 2.90)
tmin <- c(22.09, 18.65, 21.44)
tmax <- c(33.05, 30.70, 28.74)
rsds <- c(196.90, 178.50, 186.60)
# prec <- c(0.79, 0.79, 0.00, 2.90)
# tmin <- c(22.09, 22.09, 18.65, 21.44)
# tmax <- c(33.05, 33.05, 30.70, 28.74)
# rsds <- c(196.90, 196.90, 178.50, 186.60)

ncellsTotal <- 1

elevation <- 434
lat <- array(-8.25, dim = ncellsTotal)
lon <- array(-39.25, dim = ncellsTotal)
# lat <- c(  0.75,  -8.25,   0.75,  -8.25)
# lon <- c(-55.75, -55.75, -39.25, -39.25)

## Generate test forcing
rsdsAll <- tminAll <- tmaxAll <- precAll <- array(NA, dim = c(metGen$derived$nrec_in, ncellsTotal))
lonAll <- latAll <- elevationAll <- array(NA, dim = c(ncellsTotal))
for (iTime in 1:metGen$derived$nrec_in) {
  for (iCell in 1:ncellsTotal) {
    precAll[iTime, iCell] <- prec[iTime]
    tminAll[iTime, iCell] <- tmin[iTime]
    tmaxAll[iTime, iCell] <- tmax[iTime]
    rsdsAll[iTime, iCell] <- rsds[iTime]
    elevationAll[iCell] <- elevation
    lonAll[iCell] <- lon[iCell]
    latAll[iCell] <- lat[iCell]
  }
}

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
mgsetOutVars(c( "shortwave"))


## LOAD MASK/ELEVATION
elevation <- ncLoad(file = metGen$internal$ncFileNameElevation,
                    varName = metGen$settings$elevation$ncname,
                    lonlatbox = metGen$settings$lonlatbox)
mask<-elevation
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
if (!is.null(solar_geom) && dim(solar_geom$lats) == dim(mask$xyCoords$y) && dim(solar_geom$lons) == dim(mask$xyCoords$x) && dim(solar_geom$elev) == dim(elevation$Data)) {
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
  solar_geom$elevation <- elevation$Data
  
  for (ilat in 1:length(mask$xyCoords$y)) {
    print(ilat)
    solar_geom_temp <- solar_geom(lat = mask$xyCoords$y[ilat])
    solar_geom$tiny_rad_fract[,,ilat] <- solar_geom_temp$tiny_rad_fract
    solar_geom$daylength[,ilat] <- solar_geom_temp$daylength
    solar_geom$flat_potrad[,ilat] <- solar_geom_temp$flat_potrad
    for (ilon in 1:length(mask$xyCoords$x)) {
      if (!is.na(elevation$Data[ilon,ilat])) {
        solar_geom$tt_max[,ilon,ilat] <- solar_geom_tt_max0(lat = mask$xyCoords$y[ilat], elev = elevation$Data[ilon,ilat], lr = lapse_rate)
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
  # printf("Day: %d\n", iday)
  
  # ## Set Forcing for current timestep
  # precCells <- precAll[iday,]
  # tminCells <- tminAll[iday,]
  # tmaxCells <- tmaxAll[iday,]
  # rsdsCells <- rsdsAll[iday,]
  
  ## LOAD WHOLE DOMAIN FROM NETCDF
  profile$start.time.read <- Sys.time()
  inData <- readAllForcing(mask, iday)
  profile$end.time.read <- Sys.time()
  
  for (ilon in 1:length(mask$xyCoords$x)) {
    for (ilat in 1:length(mask$xyCoords$y)) {
      # lon <- lonAll[icell]
      # lat <- latAll[icell]
      # prec <- precCells[icell]
      # tmin <- tminCells[icell]
      # tmax <- tmaxCells[icell]
      # rsds <- rsdsCells[icell]
      # elevation <- elevationAll[icell]
      
      lon <- mask$xyCoords$x[ilon]
      lat <- mask$xyCoords$y[ilat]
      prec <- inData$pr[ilon, ilat,1]
      tmin <- inData$tasmin[ilon, ilat,1]
      tmax <- inData$tasmax[ilon, ilat,1]
      rsds <- inData$shortwave[ilon, ilat,1]
      
      elevation <- mask$Data[ilon, ilat]
      printf("Day:  %d, elevation %d, lon %5.2f, lat %5.2f\n", iday, elevation, lon, lat)
      
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
      yday<-metGen$derived$inYDays[iday]
      mt$ttmax0<-solar_geom$tt_max[yday,ilon,ilat]
      mt$flat_potrad<-solar_geom$flat_potrad[yday,ilat]
      mt$slope_potrad<-solar_geom$flat_potrad[yday,ilat]
      mt$daylength<-solar_geom$daylength[yday,ilat]
      tiny_radfract <- solar_geom$tiny_rad_fract[,,ilat]
      
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
      outData<-NULL
      outData$shortwave<-array(NA, dim = metGen$derived$nOutStepDay)
      # // Transfer hourly_rad to atmos structure
      for(rec in 1:metGen$derived$nOutStepDay) {
        sum = 0;
        hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
        if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
        outData$shortwave[rec] <- 0;
        for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
          if (idx_tmp <= 24) {
            idx<-idx_tmp
            outData$shortwave[rec] <- outData$shortwave[rec] + hourly_rad_prev[idx];
            # printf("rec %d, hour_local_out %d, prev_day index: %d\n", rec, hour_local_out, idx)
          } else {
            idx<-idx_tmp-24
            outData$shortwave[rec] <- outData$shortwave[rec] + hourly_rad[idx];
            # printf("rec %d, hour_local_out %d, curr_day index: %d\n", rec, hour_local_out, idx)
          }
        }
        outData$shortwave[rec] <- outData$shortwave[rec] / metGen$derived$outDt
        sum <- sum + outData$shortwave[rec]
      }
      
      ## Move data to previous timestep
      hourly_rad_prev <- hourly_rad 
      
      print(outData$shortwave)
    }
  }
}

profile$end.time.total <- Sys.time()
cat(sprintf("  Total: %.1f seconds\n", as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")))

ncellsTotal <- dim(mask$Data)[1] * dim(mask$Data)[2]
totTime<-as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")
cat(sprintf("  Total: %.1f seconds for %d x %d = %d cells. So: 100 years (67420 cels) will take: %.1f days\n", totTime, metGen$derived$nday, ncellsTotal, 
            (metGen$derived$nday * ncellsTotal), ((totTime * 67420 * 365.25 * 100) / (metGen$derived$nday * ncellsTotal) / 86400
)))

