## ISSUES
## if eg prev_pr is used than the map is moved one gridcell upwards 
## vapor pressure is too high???


rm(list=ls(all=TRUE))
# setwd("~/Documents/WORKDIRS/RProj/metGeneratoR")
shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}
library(metGeneratoR)
source('./R/Cpp_met.R')
profile<-NULL

# mgsetLonlatbox(c(92.25, 110.25, 7.25, 36.25))
# mgsetLonlatbox(c(92.25, 92.25, -8.25, -8.25))
# mgsetLonlatbox(c(92.25, 92.75, 34.25, 36.75))
mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))
mgsetPeriod(startdate = "1950-6-01", enddate = "1950-6-01")
# mgsetPeriod(startdate = "1964-12-31", enddate = "1965-01-3")
# mgsetPeriod(startdate = "1965-01-01", enddate = "1965-01-2")
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
  pr         = list(ncname = "pr",      filename = "../example_data4mtclim/Global/pr_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  tasmin     = list(ncname = "tasmin",  filename = "../example_data4mtclim/Global/tasmin_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  tasmax     = list(ncname = "tasmax",  filename = "../example_data4mtclim/Global/tasmax_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  # pressure   = list(ncname = "ps",      filename = "../example_data4mtclim/Global/ps_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  # relhum     = list(ncname = "hurs",    filename = "../example_data4mtclim/Global/hurs_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  shortwave  = list(ncname = "rsds",    filename = "../example_data4mtclim/Global/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  # longwave   = list(ncname = "rlds",    filename = "../example_data4mtclim/Global/rlds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  wind       = list(ncname = "sfcWind", filename = "../example_data4mtclim/Global/wind_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc")
))

## Define elevation file
mgsetElevation(ncname = "elevation", filename = metGen$internal$ncFileNameElevation)

## Define output variables
# mgsetOutVars(c("pr", "tas"))
# mgsetOutVars(c( "shortwave", "longwave", "tas", "pr", "pressure", "wind", "vp"))
mgsetOutVars(c( "shortwave"))


## Load elevation
elev <- ncLoad(file = metGen$internal$ncFileNameElevation,
               varName = metGen$settings$elevation$ncname,
               lonlatbox = metGen$settings$lonlatbox)

## Load mask (or base it on elevation)
mask<-elev
mask$Data[!is.na(mask$Data)]<- 1

nx <- length(mask$xyCoords$x)
ny <- length(mask$xyCoords$y)

# ## Calculate solar GEOMs as preprocessing step
lapse_rate <- 0.0065

# load("./hpc/outRadFractions_2880.Rdata")
load("./hpc/outDaylength_2880.Rdata")
load("./hpc/outFlat_potrad_2880.Rdata")
load("./hpc/outTt_max0.Rdata")

solar_geom<-NULL
# tiny_rad_fract <- aperm(outRadFractions, c(2,3,1))
# solar_geom$tiny_rad_fract <- aperm(outRadFractions, c(2,3,1))
solar_geom$daylength <- aperm(outDaylength, c(2,1))
solar_geom$flat_potrad <- aperm(outFlat_potrad, c(2,1))
solar_geom$tt_max <-  aperm(outTt_max0, c(3,1,2))
solar_geom$lats <- elev$xyCoords$y
solar_geom$lons <- elev$xyCoords$x
solar_geom$elevation <- elev$Data
# rm(outRadFractions, outDaylength, outFlat_potrad, outTt_max0)
rm(outDaylength, outFlat_potrad, outTt_max0)

# ### TODO: check if start yday is correct!!!

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
  metGen$current$timestep <- iday
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
    lon      <- mask$xyCoords$x[ilon]

    ## Calculate offset longitude
    nrOffsetSteps <- 24
    # nrOffsetSteps <- 720
    
    hour_offset <- hour_offset_int <- ceiling(ilon * (nrOffsetSteps/720))    # hour_offset<-0
  
    for (ilat in 1:length(mask$xyCoords$y)) {
      ## Select mt data for current day in year
      # mt$ttmax0       <- solar_geom$tt_max[yday, ilon, ilat]
      flat_potrad  <- solar_geom$flat_potrad[yday, ilat]
      slope_potrad <- solar_geom$flat_potrad[yday, ilat]
      daylength    <- solar_geom$daylength[yday, ilat]
      
      elevation <- mask$Data[ilon, ilat]
      if (!is.na(elevation) && !is.na(inData$pr[ilon, ilat,1])) {

        lat      <- mask$xyCoords$y[ilat]
        prec     <- inData$pr[ilon, ilat,1]
        pressure <- inData$pressure[ilon, ilat,1]
        relhum   <- inData$relhum[ilon, ilat,1] / 100 # convert to fraction
        wind <- inData$wind[ilon, ilat,1]
        tmin <- inData$tasmin[ilon, ilat,1] - 273.15
        tmax <- inData$tasmax[ilon, ilat,1] - 273.15
        shortwave <- inData$shortwave[ilon, ilat,1]
        longwave <- inData$longwave[ilon, ilat,1]
        
        if(!is.null(metGen$settings$inVar$shortwave) && !is.null(outData$shortwave)) {
          hourly_rad <- solar_geom_c(lat, yday)
          for (j in 1:nrOffsetSteps) {
            hourly_rad[j] <-  hourly_rad[j] * (shortwave*nrOffsetSteps)
          }
          
          ## Fill the arrays
          ## Copy data to previous timestep for first timestep
          if (iday == 1) {
            hourly_rad_prev <- hourly_rad 
          }

          # print("start")
          for(rec in 1:metGen$derived$nOutStepDay) {
            sum = 0;
            outData$shortwave[ilon, ilat, rec] <- 0;
            arr<-1:nrOffsetSteps
            arr_new<-shifter(arr,hour_offset_int)
            arr_new2<-arr_new[((rec-1) * metGen$derived$outDt):(rec * metGen$derived$outDt)]
            for (idx in arr_new2) {
                outData$shortwave[ilon, ilat, rec] <- outData$shortwave[ilon, ilat, rec] + hourly_rad[idx];
            }
            outData$shortwave[ilon, ilat, rec] <- outData$shortwave[ilon, ilat, rec] / metGen$derived$outDt
            sum <- sum + outData$shortwave[ilon, ilat, rec]
          }
        }
        
        ## Temperature!!
        if(!is.null(metGen$settings$inVar$tasmin) && !is.null(metGen$settings$inVar$tasmin) && !is.null(outData$tas)) {
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
              } else {
                idx<-idx_tmp-24
                outData$tas[ilon, ilat, rec] <- outData$tas[ilon, ilat, rec] + hourly_tair[idx];
              }
            }
            outData$tas[ilon, ilat, rec] <- outData$tas[ilon, ilat, rec] / metGen$derived$outDt
            sum <- sum + outData$tas[ilon, ilat, rec]
          }
        }
        # /*************************************************
        #   Precipitation
        # *************************************************/
        if(!is.null(metGen$settings$inVar$pr) && !is.null(outData$pr)) {
          if (iday == 1) {
            prec_prev <- prec
          }
          
          if(is.null(metGen$settings$inVar$pr)) {
            cat("niet gegeven!\n")
          } else {
            # cat("gegeven!\n")
            for(rec in 1:metGen$derived$nOutStepDay) {
              hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
              if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
              outData$pr[ilon, ilat, rec] <- 0;
              for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
                if (idx_tmp <= 24) {
                  outData$pr[ilon, ilat, rec] <- prec_prev;
                } else {
                  outData$pr[ilon, ilat, rec] <- prec;
                }
              }
            }
          }
        }
        
        # /**************************************
        #   Estimate Atmospheric Pressure (Pa) 
        # **************************************/
        if(!is.null(metGen$settings$inVar$pressure) && !is.null(outData$pressure)) {
          if (iday == 1) {
            pressure_prev <- pressure
          }
          
          if(is.null(metGen$settings$inVar$pressure)) {
            # cat("niet gegeven!\n")
          } else {
            # cat("gegeven!\n")
            for(rec in 1:metGen$derived$nOutStepDay) {
              hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
              if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
              outData$pressure[ilon, ilat, rec] <- 0;
              for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
                if (idx_tmp <= 24) {
                  outData$pressure[ilon, ilat, rec] <- pressure_prev;
                } else {
                  outData$pressure[ilon, ilat, rec] <- pressure;
                }
              }
            }
          }
        }
        
        # /*************************************************
        #   Wind Speed
        # *************************************************/
        if(!is.null(metGen$settings$inVar$wind) && !is.null(outData$wind)) {
          if (iday == 1) {
            wind_prev <- wind
          }
          
          if(is.null(metGen$settings$inVar$wind)) {
            # cat("niet gegeven!\n")
          } else {
            # cat("gegeven!\n")
            for(rec in 1:metGen$derived$nOutStepDay) {
              hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
              if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
              outData$wind[ilon, ilat, rec] <- 0;
              for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
                if (idx_tmp <= 24) {
                  outData$wind[ilon, ilat, rec] <- wind_prev;
                } else {
                  outData$wind[ilon, ilat, rec] <- wind;
                }
              }
            }
          }
        }
        
        # /*************************************************
        #   Longwave
        # *************************************************/
        if(!is.null(metGen$settings$inVar$pressure) && !is.null(outData$pressure)) {
          if (iday == 1) {
            longwave_prev <- longwave
          }
          
          if(is.null(metGen$settings$inVar$longwave)) {
            # cat("niet gegeven!\n")
          } else {
            # cat("gegeven!\n")
            for(rec in 1:metGen$derived$nOutStepDay) {
              hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
              if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
              outData$longwave[ilon, ilat, rec] <- 0;
              for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
                if (idx_tmp <= 24) {
                  outData$longwave[ilon, ilat, rec] <- longwave_prev;
                } else {
                  outData$longwave[ilon, ilat, rec] <- longwave;
                }
              }
            }
          }
        }
        
        # /*************************************************
        #   Vapor pressure
        # *************************************************/
        if(!is.null(metGen$settings$inVar$relhum) && !is.null(outData$vp)) {
          
          if (iday == 1) {
            relhum_prev <- relhum
          }
          
          if(!is.null(metGen$settings$inVar$relhum) && !is.null(metGen$settings$outVars$tas)) {
            for(rec in 1:metGen$derived$nOutStepDay) {
              hour_local_out <- (rec-1)*metGen$derived$outDt - hour_offset_int;
              if ((0 - hour_offset_int) < 0) hour_local_out <- hour_local_out + 24;
              outData$vp[ilon, ilat, rec] <- 0;
              for (idx_tmp in (hour_local_out+1):(((hour_local_out-1)+metGen$derived$outDt)+1)) {
                if (idx_tmp <= 24) {
                  outData$vp[ilon, ilat, rec] <- relhum_prev * svp(outData$tas[ilon, ilat, rec]) / 100
                } else {
                  outData$vp[ilon, ilat, rec] <- relhum * svp(outData$tas[ilon, ilat, rec]) / 100
                }
              }
            }
          }
        }
        
        ## Move data to previous timestep
        hourly_rad_prev <- hourly_rad 
        # hourly_tair_prev <- hourly_tair 
        prec_prev <- prec 
        # pressure_prev <- pressure
        # wind_prev <- wind
        # longwave_prev <- longwave
        # relhum_prev <- relhum
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
