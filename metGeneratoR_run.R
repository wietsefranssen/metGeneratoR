## ISSUES/TODO
## if eg prev_pr is used than the map is moved one gridcell upwards 
## vapor pressure is too high???
## swdown can only be done globally: fix this!
## check if start yday is correct
## check base_offset it is now right for 6hourly calculations but maybe not for others like 3 hourly
## unitconversion

## Cleanup
rm(list=ls(all=TRUE))

## Load the library
library(metGeneratoR)

metGen$settings$timeUnitsName <- "unit"
metGen$settings$uselessDim <- T

## Override standard settings
#mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))
mgsetLonlatbox(c(-83.25,-32.75,-57.75,14.25))

# mgsetPeriod(startdate = "2010-1-1", enddate = "2010-1-2")
mgsetPeriod(startdate = "1993-2-1", enddate = "1993-2-2")
# mgsetLonlatbox(c(92.25, 110.25, 7.25, 36.25))
# 
# ## Set the period to run
# mgsetPeriod(startdate = "1998-1-1", enddate = "1998-1-31")

mgsetInVars(list(
  swdown         = list(ncname = "SWdown",      filename = "~/input/199302/E01/E01_SWdown.nc"),
  lwdown         = list(ncname = "LWdown",      filename = "~/input/199302/E01/E01_LWdown.nc")
  # tasmin     = list(ncname = "tasmin",  filename = "../example_data4mtclim/Global/tasmin_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  # tasmax     = list(ncname = "tasmax",  filename = "../example_data4mtclim/Global/tasmax_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  # psurf      = list(ncname = "ps",      filename = "../example_data4mtclim/Global/ps_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  # # relhum     = list(ncname = "hurs",    filename = "../example_data4mtclim/Global/hurs_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  # qair       = list(ncname = "Qair",    filename = "../example_data4mtclim/Qair_1998.nc"),
  # swdown     = list(ncname = "rsds",    filename = "../example_data4mtclim/Global/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc")
  # lwdown     = list(ncname = "rlds",    filename = "../example_data4mtclim/Global/rlds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  # wind       = list(ncname = "sfcWind", filename = "../example_data4mtclim/Global/wind_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc")
))


mgsetInDt(inDt = 24)
mgsetOutDt(outDt = 6)

mgsetOutVars(c( "swdown","lwdown"))

mgsetOutName("./output/<VAR>_6hourly2_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_<SYEAR>_local.nc")

## Run metGen
metGenRun()

