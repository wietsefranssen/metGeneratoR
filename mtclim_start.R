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
mt <- mtclim_init(have_dewpt, have_shortwave, elevation, slope, aspect, ehoriz, whoriz, annual_prcp,
                  lat, dmy, prec, tmax, tmin, vp, hourlyrad,
                  tiny_radfract,
                  p, mtclim_data)

mt<-calc_tair(mt)

mt<-calc_prcp(mt)

mt<-snowpack(mt)

profile$start.time.run <- Sys.time()
mt1<-calc_srad_humidity_iterative(mt, options)
profile$end.time.run <- Sys.time()
cat(sprintf("  Times (run): %.1f seconds\n",
            as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))

profile$start.time.run <- Sys.time()
mt2<-calc_tiny_radfract(mt, options)
profile$end.time.run <- Sys.time()
cat(sprintf("  Times (run): %.1f seconds\n",
            as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))

save(mt2, file = "mt2.Rdata")












calc_rest_1t <- function(mt, options)
{
  p <- mt$p
  ctrl <- mt$ctrl
  data <- mt$mtclim_data
  # tiny_radfract <- mt$tiny_radfract
  ttmax0 <- mt$ttmax0
  flat_potrad <- mt$flat_potrad
  slope_potrad <- mt$slope_potrad
  daylength <- mt$daylength
  
  constants <- metGen$constants
  
  # /* number of simulation days */
  ndays <- 1
  
  # /* local array memory allocation */
  # dtr <- array(0, dim = c(ndays))
  # parray <- array(0, dim = c(ndays))
  # t_fmax <- array(0, dim = c(ndays))
  # tdew <- array(0, dim = c(ndays))
  # pet <- array(0, dim = c(ndays))
  # pva <- array(0, dim = c(ndays))
  # tdew_save <- array(0, dim = c(ndays))
  # pva_save <- array(0, dim = c(ndays))
  window <- array(0, dim = c(ndays+90))
  
  # /* calculate diurnal temperature range for transmittance calculations */
  # for (i in 1:ndays) {
  tmax = data$tmax;
  tmin = data$tmin;
  if (tmax < tmin) 
    tmax = tmin;
  dtr = tmax-tmin;
  # }
  
  # /* smooth dtr array: After Bristow and Campbell, 1984 */
  sm_dtr <- dtr
  
  # /* calculate the annual total precip */
  ann_prcp <- data$s_prcp * 365.25;
  if (ann_prcp == 0.0) ann_prcp <- 1.0;
  effann_prcp = data$s_prcp * 365.25;
  if (effann_prcp < 8.0) effann_prcp = 8.0;
  parray <- effann_prcp;
  
  # /* b parameter, and t_fmax not varying with Tdew, so these can be
  # calculated once, outside the iteration between radiation and humidity
  # estimates. Requires storing t_fmax in an array. */
  # /* b parameter from 30-day average of DTR */
  b = constants$B0 + constants$B1 * exp(-constants$B2 * sm_dtr);
  
  # /* proportion of daily maximum transmittance */
  t_fmax = 1.0 - 0.9 * exp(-b * dtr^constants$C);
  
  # /* correct for precipitation if this is a rain day */
  if (data$prcp > options$SW_PREC_THRESH) t_fmax <- t_fmax * constants$RAIN_SCALAR;
  data$s_tfmax <- t_fmax;
  
  # /* Initial values of vapor pressure, etc */
  if (ctrl$indewpt) {
    # /* Observed Tdew supplied */
    tdew <- data$tdew;
  } else {
    # /* Estimate Tdew */
    tdew <- data$s_tmin
  }
  if (ctrl$invp) {
    # /* Observed vapor pressure supplied */
    pva = data$s_hum;
  } else {
    # /* convert dewpoint to vapor pressure */
    pva = svp(tdew);
  }
  
  # /* Other values needed for srad_humidity calculation */
  pa = atm_pres(p$site_elev);
    yday = data$yday
    data$s_dayl = daylength
    tdew_save = tdew
    pva_save = pva

  # /* Initial estimates of solar radiation, cloud fraction, etc. */
  initVals <- compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, 
                                            slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options);
  pet <- initVals$pet
  tdew <- initVals$tdew
  pva <- initVals$pva
  data <- initVals$data
  # # /* estimate annual PET */
  sum_pet = 0.0;
    sum_pet <- sum_pet + pet
  ann_pet = (sum_pet/ndays) * 365.25;
  
  VP_ITER_ALWAYS <- 1
  VP_ITER_ANNUAL <- 2
  VP_ITER_CONVERGE <- 3
  
  # # /* Reset humidity terms if no iteration desired */
  if (ctrl$indewpt || ctrl$invp || (options$VP_ITER == VP_ITER_ANNUAL && ann_pet/ann_prcp >= 2.5) ) {
      tdew = tdew_save
      pva = pva_save
  }
  
  # /* Set up srad-humidity iterations */
  if (options$VP_ITER == VP_ITER_ALWAYS || (options$VP_ITER == VP_ITER_ANNUAL && ann_pet/ann_prcp >= 2.5) || options$VP_ITER == VP_ITER_CONVERGE) {
    # //printf("Using arid-climate humidity algorithm\n");
    if (options$VP_ITER == VP_ITER_CONVERGE) {
      max_iter = 100;
    } else {
      max_iter = 2;
    }
  } else {
    # //printf("Using Tdew=Tmin humidity algorithm\n");
    max_iter = 1;
  }
  
  # /* srad-humidity iterations */
  tol = 0.01;
  iter = 1;
  rmse_tdew = tol+1;
  while (rmse_tdew > tol && iter < max_iter) {
    # print(iter)
    update_pva = 1;
    for (i in 1:ndays) {
      tdew_save[i] = tdew[i];
    }
    initVals<-compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, 
                                            slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options);
    pet <- initVals$pet
    tdew <- initVals$tdew
    pva <- initVals$pva
    data <- initVals$data
    rmse_tdew = 0;
    for (i in 1:ndays) {
      rmse_tdew<-rmse_tdew + ((tdew[i]-tdew_save[i])*(tdew[i]-tdew_save[i]))
    }
    rmse_tdew<-rmse_tdew / ndays;
    rmse_tdew = rmse_tdew^0.5;
    iter <- iter + 1
  }
  
  # /* save humidity in output data structure */
  if (ctrl$outhum) {
    if (!ctrl$invp) {
      for (i in 1:ndays) {
        data$s_hum[i] = pva[i];
      }
    }
  } else {
    # /* output humidity as vapor pressure deficit (Pa) */
    for (i in 1:ndays) {
      # /* calculate saturated VP at tday */
      # /* start vic_change */
      # /* pvs = 610.7 * exp(17.38 * data$s_tday[i]/(239.0+data$s_tday[i])); */
      pvs = svp(data$s_tday[i]);
      # /* end vic_change */
      vpd = pvs - pva[i];
      if (vpd < 0.0) vpd = 0.0;
      data$s_hum[i] = vpd;
    }
  }
  
  mt$p <- p
  mt$ctrl <- ctrl
  mt$mtclim_data <- data
  # mt$tiny_radfract <- tiny_radfract
  mt$ttmax0 <- ttmax0
  mt$flat_potrad <-flat_potrad
  mt$slope_potrad <- slope_potrad
  mt$daylength <- daylength
  
  return(mt)
} # /* end of calc_srad_humidity_iterative() */














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

profile$start.time.run <- Sys.time()
for (i in 1:(360*5)) mt_t1<-calc_rest_1t(mt2, options)
profile$end.time.run <- Sys.time()
cat(sprintf("  New Times (run): %.1f seconds\n",
            as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs")))
