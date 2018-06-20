rm(list=ls(all=TRUE))

library(metGeneratoR)
profile<-NULL
source("./R/temp_function.R")
# mgsetPeriod(startdate = "1950-01-03", enddate = "1950-01-05")
# mgsetPeriod(startdate = "1964-12-31", enddate = "1965-01-03")
mgsetPeriod(startdate = "1965-01-01", enddate = "1965-01-03")
mgsetInDt(24) # Set N hours per timestep
mgsetOutDt(6) # Set N hours per timestep

have_dewpt <- F
elevation <- 434
annual_prcp <- 583.219500000
lat <- -8.25
lon <- -39.25
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
param_set_TYPE_SHORTWAVE_SUPPLIED <- T
options <- NULL
options$SW_PREC_THRESH <- 0
options$MTCLIM_SWE_CORR <- 0
options$LW_CLOUD <- 1
options$VP_ITER <- 1

outData<-NULL

theta_l = -30
theta_s = lon;
theta_l = 0.25
theta_s = 0.25;
hour_offset = (theta_l-theta_s)*24/360;
if (hour_offset < 0) {
  hour_offset_int = floor(hour_offset-0.5);
  cat("sd" )
} else {
  hour_offset_int = floor(hour_offset+0.5);
}
hour_offset <- hour_offset - hour_offset_int #// hour_offset is now the distance from the center of local time zone
print(hour_offset_int)






metGen$constants<-setConstants()
constants <- metGen$constants

# /*************************************************
#   Shortwave, part 1.
# *************************************************/
hourlyrad <- array(NA, dim = c(metGen$derived$nday*24))
if (param_set_TYPE_SHORTWAVE_SUPPLIED) {
  have_shortwave = 1; # // flag for MTCLIM
  for (day in 1:metGen$derived$nday) {
    for (hour in 1:24) {
      if(metGen$derived$inDt == 24) {
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
stop()
###################### START MTCLIM WRAPPER

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
for (i in 1:(420)) mt_t2<-calc_rest_1t(mt2, options)
# for (i in 1:(67420)) mt_t2<-calc_rest_1t(mt2, options)
profile$end.time.run <- Sys.time()
cat(sprintf("  New Times (run): %.1f seconds\n",
            as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))

ctrl<-mt$ctrl
mtclim_data<-mt_t2
yday
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

rec<-1




outData$shortwave<-array(NA, dim = metGen$derived$outDt)
# // Transfer hourlyrad to atmos structure
for(rec in 1:metGen$derived$nrec_out) {
  sum = 0;
  hour = (rec-1)*metGen$derived$outDt - hour_offset_int;
  if ((0 - hour_offset_int) < 0) hour <- hour + 24;
  outData$shortwave[rec] <- 0;
  for (idx in (hour+1):(((hour-1)+metGen$derived$outDt)+1)) {
    outData$shortwave[rec] <- outData$shortwave[rec] + hourlyrad[idx];
    printf("rec %d, hour %d, start hour: %d (corresponds with index: %d)\n", rec,hour, idx-24-1, idx)
  }
  outData$shortwave[rec] <- outData$shortwave[rec] / metGen$derived$outDt
  sum <- sum + outData$shortwave[rec]
}
