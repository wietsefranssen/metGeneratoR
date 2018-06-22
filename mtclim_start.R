rm(list=ls(all=TRUE))

library(metGeneratoR)
profile<-NULL
source("./R/temp_function.R")
# mgsetPeriod(startdate = "1950-01-03", enddate = "1950-01-05")
# mgsetPeriod(startdate = "1964-12-31", enddate = "1965-01-03")
mgsetPeriod(startdate = "1965-01-01", enddate = "1965-01-03")
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


ncellsTotal <- 5
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

## Calculate solar GEOMs as preprocessing step
lapse_rate <- 0.0065
solar_geom_file <- "./solar_geom.Rdata"
if(file.exists(solar_geom_file)) {
  load(file = solar_geom_file)
  ## Check the file:
  if (length(solar_geom$lat) == length(lat) &&
      solar_geom$lat == lat && solar_geom$elev == elevation && solar_geom$lapse_rate == lapse_rate) {
    print("Good!") 
  } else {
    print("Error in file! Regenerating solar_geom file...") 
    solar_geom<-solar_geom(lat = lat)
    solar_geom$tt_max<-solar_geom_tt_max0(lat = lat, elev = elevation, lr = lapse_rate)
    solar_geom$lat <- lat
    solar_geom$elev <- elevation
    solar_geom$lapse_rate <- lapse_rate
    save(solar_geom, file = solar_geom_file)
  }
}

### THE MAIN LOOP
profile$start.time.total <- Sys.time()
for (iday in 1:metGen$derived$nday) {
  metGen$current$timestep <- 1
  # printf("Day: %d\n", iday)
  
  ## Set Forcing for current timestep
  precCells <- precAll[iday,]
  tminCells <- tminAll[iday,]
  tmaxCells <- tmaxAll[iday,]
  rsdsCells <- rsdsAll[iday,]
  
  for (icell in 1:ncellsTotal) {
    lon <- lonAll[icell]
    lat <- latAll[icell]
    prec <- precCells[icell]
    tmin <- tminCells[icell]
    tmax <- tmaxCells[icell]
    rsds <- rsdsCells[icell]
    
    elevation <- elevationAll[icell]
    printf("Day:  %d, Cell: %d, elevation %d, lon %5.2f, lat %5.2f\n", iday, icell, elevation, lon, lat)
    
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
    hourlyrad <- array(NA, dim = c(24))
    if (param_set_TYPE_SHORTWAVE_SUPPLIED) {
      if(metGen$derived$inDt == 24) {
        hourlyrad <- rep(rsds, 24)
      } else {
        for (hour in 1:24) {
          hourlyrad[hour] = rsds[hour];
        }
      }
    }
    ###################### START MTCLIM WRAPPER
    printf("doing mtclim_init...\n")
    mt <- mtclim_init(have_dewpt, param_set_TYPE_SHORTWAVE_SUPPLIED, elevation, 0,0,0,
                      lat, prec, tmax, tmin, vp, metGen$derived$inYDays[iday], hourlyrad,
                      tiny_radfract,
                      p, mtclim_data)
    
    mt<-calc_tair(mt)
    
    mt<-calc_prcp(mt)
    
    mt<-snowpack(mt)
    
    ## Select mt data for current day in year
    yday<-metGen$derived$inYDays[iday]
    mt$ttmax0<-solar_geom$tt_max[yday]
    mt$flat_potrad<-solar_geom$flat_potrad[yday]
    mt$slope_potrad<-solar_geom$flat_potrad[yday]
    mt$daylength<-solar_geom$daylength[yday]
    tiny_radfract <- solar_geom$tiny_rad_fract
    
    profile$start.time.run <- Sys.time()
    for (i in 1:(1)) mt_t2<-calc_rest_1t(mt, options)
    # for (i in 1:(67420)) mt_t2<-calc_rest_1t(mt2, options)
    profile$end.time.run <- Sys.time()
    cat(sprintf("  Times (run): %.1f seconds\n",
                as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))
    ctrl<-mt$ctrl
    mtclim_data<-mt_t2
    hourly_rad <- mtclim_to_vic(hour_offset, yday, tiny_radfract, mt$ctrl, 
                                mt_t2, tskc, vp, fdir)
    
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
          hourlyrad[(day-1)*24+hour] = local_forcing_data[SHORTWAVE][(day-1)*24+hour];
        }
      }
    }
    
    ## Copy data to previous timestep for first timestep
    if (iday == 1) {
      hourlyrad_prev <- hourlyrad 
    }
    
    
    ## Fill the arrays
    outData<-NULL
    outData$shortwave<-array(NA, dim = metGen$derived$nOutStepDay)
    # // Transfer hourlyrad to atmos structure
    for(rec in 1:metGen$derived$nOutStepDay) {
      sum = 0;
      hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
      if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
      outData$shortwave[rec] <- 0;
      for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
        if (idx_tmp <= 24) {
          idx<-idx_tmp
          outData$shortwave[rec] <- outData$shortwave[rec] + hourlyrad_prev[idx];
          # printf("rec %d, hour_local_out %d, prev_day index: %d\n", rec, hour_local_out, idx)
        } else {
          idx<-idx_tmp-24
          outData$shortwave[rec] <- outData$shortwave[rec] + hourlyrad[idx];
          # printf("rec %d, hour_local_out %d, curr_day index: %d\n", rec, hour_local_out, idx)
        }
      }
      outData$shortwave[rec] <- outData$shortwave[rec] / metGen$derived$outDt
      sum <- sum + outData$shortwave[rec]
    }
    # print(outData$shortwave)
    
    ## Move data to previous timestep
    hourlyrad_prev <- hourlyrad 
  }
}

profile$end.time.total <- Sys.time()
cat(sprintf("  Total: %.1f seconds\n", as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")))

totTime<-as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")
cat(sprintf("  Total: %.1f seconds for %d x %d = %d cells. So: 100 yeas (67420 cels) will take: %.1f days\n", totTime, metGen$derived$nday ,ncellsTotal, (metGen$derived$nday * ncellsTotal), (
  (totTime * 67420 * 365.25 * 100) / (metGen$derived$nday * ncellsTotal) / 86400
)))

