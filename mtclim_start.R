rm(list=ls(all=TRUE))

library(metGeneratoR)
profile<-NULL
source("./R/temp_function.R")
# mgsetPeriod(startdate = "1950-01-03", enddate = "1950-01-05")
# mgsetPeriod(startdate = "1964-12-31", enddate = "1965-01-03")
mgsetPeriod(startdate = "1965-01-01", enddate = "1965-01-01")
mgsetNHourPerStep(6) # Set N hours per timestep
metGen$settings$nHourPerStepIn <- 24

have_dewpt <- F
elevation <- 434
annual_prcp <- 583.219500000
lat <- -8.25
# lat <- 79.25
slope <- 0
aspect <- 0
ehoriz <- whoriz <- 0

prec <- c(0.79, 0.00, 2.90)
tmin <- c(22.09, 18.65, 21.44)
tmax <- c(33.05, 30.70, 28.74)
rsds <- c(196.90, 178.50, 186.60)
prec <- c(0.79, 0.79, 0.00, 2.90)
tmin <- c(22.09, 22.09, 18.65, 21.44)
tmax <- c(33.05, 33.05, 30.70, 28.74)
rsds <- c(196.90, 196.90, 178.50, 186.60)
# prec <- c(0.79, 0.00)
# tmin <- c(22.09, 18.65)
# tmax <- c(33.05, 30.70)
# rsds <- c(196.90, 178.50)
param_set_TYPE_SHORTWAVE_SUPPLIED <- T
options <- NULL
options$SW_PREC_THRESH <- 0
options$MTCLIM_SWE_CORR <- 0
options$LW_CLOUD <- 1
options$VP_ITER <- 1

metGen$constants<-setConstants()
constants <- metGen$constants

# /*************************************************
#   Shortwave, part 1.
# *************************************************/
hourlyrad <- array(NA, dim = metGen$derived$nday)
if (param_set_TYPE_SHORTWAVE_SUPPLIED) {
  have_shortwave = 1; # // flag for MTCLIM
  for (day in 1:metGen$derived$nday) {
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
# 
## Run mtclim_init
# mt <- mtclim_init(have_dewpt, have_shortwave, elevation, slope, aspect, ehoriz, whoriz, annual_prcp,
#                   lat, dmy, prec, tmax, tmin, vp, hourlyrad,
#                   tiny_radfract,
#                   p, mtclim_data)
# 
# mt<-calc_tair(mt)
# 
# mt<-calc_prcp(mt)
# 
# mt<-snowpack(mt)
# 
# profile$start.time.run <- Sys.time()
# mt1<-calc_srad_humidity_iterative(mt, options)
# profile$end.time.run <- Sys.time()
# cat(sprintf("  Times (run): %.1f seconds\n",
#             as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))
# 
# profile$start.time.run <- Sys.time()
# mt2<-calc_tiny_radfract(mt, options)
# profile$end.time.run <- Sys.time()
# cat(sprintf("  Times (run): %.1f seconds\n",
#             as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))
# 
# save(mt2, file = "mt2.Rdata")
# 
# 

load("mt2.Rdata")
tiny_radfract <- mt2$tiny_radfract
mt2$tiny_radfract<-NULL

yday<-metGen$derived$inYDays
mt<-mt2
mt$ttmax0<-mt$ttmax0[yday]
mt$flat_potrad<-mt$flat_potrad[yday]
mt$slope_potrad<-mt$slope_potrad[yday]
mt$daylength<-mt$daylength[yday]


profile$start.time.run <- Sys.time()
for (i in 1:(360*5)) mt3<-calc_rest2(mt2, options)
profile$end.time.run <- Sys.time()
cat(sprintf("  Old Times (run): %.1f seconds\n",
            as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))

mt2<-mt
profile$start.time.run <- Sys.time()
for (i in 1:(76420)) mt_t2<-calc_rest_1t(mt2, options)
profile$end.time.run <- Sys.time()
cat(sprintf("  New Times (run): %.1f seconds\n",
            as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))
