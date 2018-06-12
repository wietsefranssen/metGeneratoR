mtclim_getSettings <- function() {

  ## Init the settings datatype
  settings <- initSettings(startdate = "1950-01-01",
                           enddate = "1950-1-31",
                           outstep = 6,
                           lonlatbox = c(92.25, 110.25, 7.25, 36.25),
                           outfile = "example_output_Mekong_",
                           outperyear = TRUE)

  ## Some system settings
  settings$system$nCores <- 2
  settings$system$maxMem <- 0.0040 # in Gb

  ## Return the location of the Example NetCDF-files
  ncFileNameElevation  <- system.file("extdata", "elevation_Mekong.nc4", package = "mtclimR")
  ncFileNamePr         <- system.file("extdata", "pr_Mekong.nc4", package = "mtclimR")
  ncFileNameTasmin     <- system.file("extdata", "tasmin_Mekong.nc4", package = "mtclimR")
  ncFileNameTasmax     <- system.file("extdata", "tasmax_Mekong.nc4", package = "mtclimR")
  ncFileNameWind       <- system.file("extdata", "wind_Mekong.nc4", package = "mtclimR")

  ## Input variables
  ## Comment out the ones you dont want to include
  settings <- setInputVars(settings,list(
    pr         = list(ncFileName = ncFileNamePr,        ncName = "pr",        vicIndex = 9,   alma = FALSE),
    tasmin     = list(ncFileName = ncFileNameTasmin,    ncName = "tasmin",    vicIndex = 17),
    tasmax     = list(ncFileName = ncFileNameTasmax,    ncName = "tasmax",    vicIndex = 16),
    wind       = list(ncFileName = ncFileNameWind,      ncName = "wind",      vicIndex = 20)
  ))
  settings$elevation <- list(ncFileName = ncFileNameElevation, ncName = "elevation")

  ## Output variables
  ## Comment out the ones you dont want to include
  settings$outputVars <- list(
    pr         = list(VICName = "OUT_PREC",       units = "mm",        longName = "incoming precipitation"),
    tas        = list(VICName = "OUT_AIR_TEMP",   units = "C",         longName = "air temperature"),
    shortwave  = list(VICName = "OUT_SHORTWAVE",  units = "W m-2",     longName = "incoming shortwave"),
    longwave   = list(VICName = "OUT_LONGWAVE",   units = "W m-2",     longName = "incoming longwave"),
    pressure   = list(VICName = "OUT_PRESSURE",   units = "kPa",       longName = "near surface atmospheric pressure"),
    qair       = list(VICName = "OUT_QAIR",       units = "kg kg-1",   longName = "specific humidity"),
    vp         = list(VICName = "OUT_VP",         units = "kPa",       longName = "near surface vapor pressure"),
    rel_humid  = list(VICName = "OUT_REL_HUMID",  units = "fraction",  longName = "relative humidity"),
    density    = list(VICName = "OUT_DENSITY",    units = "kg m-3",    longName = "near-surface atmospheric density"),
    wind       = list(VICName = "OUT_WIND",       units = "m s-1",     longName = "near surface wind speed")
  )

  return(settings)
}

# ///***** Forcing Variable Types *****/
#   //#define N_FORCING_TYPES 24
#   //#define AIR_TEMP   0 /* air temperature per time step [C] (ALMA_INPUT: [K]) */
#   //#define ALBEDO     1 /* surface albedo [fraction] */
#   //#define CHANNEL_IN 2 /* incoming channel flow [m3] (ALMA_INPUT: [m3/s]) */
#   //#define CRAINF     3 /* convective rainfall [mm] (ALMA_INPUT: [mm/s]) */
#   //#define CSNOWF     4 /* convective snowfall [mm] (ALMA_INPUT: [mm/s]) */
#   //#define DENSITY    5 /* atmospheric density [kg/m3] */
#   //#define LONGWAVE   6 /* incoming longwave radiation [W/m2] */
#   //#define LSRAINF    7 /* large-scale rainfall [mm] (ALMA_INPUT: [mm/s]) */
#   //#define LSSNOWF    8 /* large-scale snowfall [mm] (ALMA_INPUT: [mm/s]) */
#   //#define PREC       9 /* total precipitation (rain and snow) [mm] (ALMA_INPUT: [mm/s]) */
#   //#define PRESSURE  10 /* atmospheric pressure [kPa] (ALMA_INPUT: [Pa]) */
#   //#define QAIR      11 /* specific humidity [kg/kg] */
#   //#define RAINF     12 /* rainfall (convective and large-scale) [mm] (ALMA_INPUT: [mm/s]) */
#   //#define REL_HUMID 13 /* relative humidity [fraction] */
#   //#define SHORTWAVE 14 /* incoming shortwave [W/m2] */
#   //#define SNOWF     15 /* snowfall (convective and large-scale) [mm] (ALMA_INPUT: [mm/s]) */
#   //#define TMAX      16 /* maximum daily temperature [C] (ALMA_INPUT: [K]) */
#   //#define TMIN      17 /* minimum daily temperature [C] (ALMA_INPUT: [K]) */
#   //#define TSKC      18 /* cloud cover fraction [fraction] */
#   //#define VP        19 /* vapor pressure [kPa] (ALMA_INPUT: [Pa]) */
#   //#define WIND      20 /* wind speed [m/s] */
#   //#define WIND_E    21 /* zonal component of wind speed [m/s] */
#   //#define WIND_N    22 /* meridional component of wind speed [m/s] */
#   //#define SKIP      23 /* place holder for unused data columns */
