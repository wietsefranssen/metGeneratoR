rm(list=ls(all=TRUE))

library(metGeneratoR)
profile<-NULL

# mgsetPeriod(startdate = "1950-01-03", enddate = "1950-01-05")
# mgsetPeriod(startdate = "1964-12-31", enddate = "1965-01-03")
mgsetPeriod(startdate = "1965-01-01", enddate = "1965-01-03")
mgsetNHourPerStep(6) # Set N hours per timestep
metGen$settings$nHourPerStepIn <- 24

have_dewpt <- F
elevation <- 434
annual_prcp <- 583.219500000
lat <- -8.25
slope <- 0
aspect <- 0
ehoriz <- whoriz <- 0

prec <- c(0.79, 0.00, 2.90)
tmin <- c(22.09, 18.65, 21.44)
tmax <- c(33.05, 30.70, 28.74)
rsds <- c(196.90, 178.50, 186.60)
param_set_TYPE_SHORTWAVE_SUPPLIED <- F
options <- NULL
options$SW_PREC_THRESH <- 0
options$MTCLIM_SWE_CORR <- 0
options$LW_CLOUD <- 1
options$VP_ITER <- 1
ndays <- length(tmin)

metGen$constants<-setConstants()
constants <- metGen$constants

# /*************************************************
#   Shortwave, part 1.
# *************************************************/
hourlyrad <- array(NA, dim = ndays)
if (param_set_TYPE_SHORTWAVE_SUPPLIED) {
  have_shortwave = 1; # // flag for MTCLIM
  for (day in 1:ndays) {
    for (hour in 1:24) {
      if(metGen$settings$nHourPerStepIn == 24) {
        hourlyrad[(day-1)*24+hour] = rsds[day];
      }
      else {
        hourlyrad[(day-1)*24+hour] = rsds[(day-1)*24+hour];
      }
    }
  }
} else {
  have_shortwave = 0;
}

## Run mtclim_init
mt <- mtclim_init(have_dewpt, have_shortwave, elevation, slope, aspect, ehoriz, whoriz, annual_prcp,
                  lat, dmy, prec, tmax, tmin, vp, hourlyrad,
                  tiny_radfract, ctrl,
                  p, mtclim_data)

mt<-calc_tair(mt)

mt<-calc_prcp(mt)

mt<-snowpack(mt)

profile$start.time.run <- Sys.time()
mt<-calc_srad_humidity_iterative(mt, options)
profile$end.time.run <- Sys.time()

cat(sprintf("  Times (run): %.1f seconds\n",
            as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))
