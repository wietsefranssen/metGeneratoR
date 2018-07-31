## ISSUES/TODO
## if eg prev_pr is used than the map is moved one gridcell upwards 
## vapor pressure is too high???
## shortwave can only be done globally: fix this!
## check if start yday is correct
## check base_offset it is now right for 6hourly calculations but maybe not for others like 3 hourly
## unitconversion

rm(list=ls(all=TRUE))

library(metGeneratoR)

profile<-NULL

# mgsetLonlatbox(c(92.25, 110.25, 7.25, 36.25))
# mgsetLonlatbox(c(92.25, 92.25, -8.25, -8.25))
# mgsetLonlatbox(c(92.25, 92.75, 34.25, 36.75))
mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))
mgsetPeriod(startdate = "1998-6-1", enddate = "1998-6-1")
# mgsetPeriod(startdate = "1965-01-01", enddate = "1965-06-2")
mgsetInDt(24) # Set N hours per timestep
mgsetOutDt(6) # Set N hours per timestep

metGen$constants<-setConstants()
constants <- metGen$constants

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
# mgsetOutVars(c("shortwave", "tas", "vp"))
mgsetOutVars(c( "shortwave", "longwave", "tas", "pr", "pressure", "wind", "vp"))

## Load elevation
elev <- ncLoad(file = metGen$internal$ncFileNameElevation,
               varName = metGen$settings$elevation$ncname,
               lonlatbox = metGen$settings$lonlatbox)

## Load mask (or base it on elevation)
mask<-elev
mask$Data[!is.na(mask$Data)]<- 1

nx <- length(mask$xyCoords$x)
ny <- length(mask$xyCoords$y)

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
  radfrac<-rad_map_final_cr(metGen$derived$nOutStepDay, yday, nx_parts = 720, gmt_float = 0)

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
      outData$longwave[, , rec] <- inData$longwave[, ,1]
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
    outData$tas<- set_max_min_lonlat_cr(inData$tasmin[,,1], inData$tasmax[,,1], yday, metGen$derived$nOutStepDay)
  }
  outData$tas<- outData$tas - 273.15
  # /*************************************************
  #   Vapor pressure
  # *************************************************/
  inData$relhum <- inData$relhum / 100
  if(!is.null(metGen$settings$inVar$relhum) && !is.null(outData$vp)) {
    for(rec in 1:metGen$derived$nOutStepDay) {
      outData$vp <- set_vp_cr(outData$tas, inData$relhum[,,1], nx, ny, metGen$derived$nOutStepDay)
    }
  }

  ## ADD OUTPUT TO NETCDF
  for (var in names(metGen$settings$outVars)) {
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

## Stats
ncellsTotal <- sum(mask$Data, na.rm = TRUE)
profile$end.time.total <- Sys.time()
totTime<-as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")
cat(sprintf("  Total: %.1f seconds for %d x %d = %d cells. So: 100 years (67420 cels) will take: %.1f days\n", totTime, metGen$derived$nday, ncellsTotal, 
            (metGen$derived$nday * ncellsTotal), ((totTime * 94742 * 365.25 * 100) / (metGen$derived$nday * ncellsTotal) / 86400
            )))
