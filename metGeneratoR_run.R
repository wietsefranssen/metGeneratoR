## Cleanup
rm(list=ls(all=TRUE))

## Load the library
library(metGeneratoR)

## Init standard settings
mgsetInit()

## Override standard settings

# mgsetLonlatbox(c(92.25, 92.75, 34.25, 36.75))
mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))

mgsetPeriod(startdate = "1998-1-1", enddate = "1998-1-2")

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

mgsetOutVars(c( "shortwave", "longwave", "tas", "pr", "pressure", "wind", "vp"))
mgsetOutVars(c( "shortwave", "tas"))

mgsetOutName("./output/<var>/<VAR>_6hourly_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_<SYEAR>.nc")

## Run metGen
metGenRun()
# metGenRunDemo()


