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
  sky_prop <- mt$sky_prop
  constants <- metGen$constants
  
  # /* local array memory allocation */
  pet <- NULL

  # /* calculate diurnal temperature range for transmittance calculations */
  tmax = data$tmax;
  tmin = data$tmin;
  if (tmax < tmin) 
    tmax = tmin;
  dtr = tmax-tmin;

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
  initVals <- compute_srad_humidity_onetime_t1(ctrl, data, tdew, pva, ttmax0, flat_potrad, 
                                            slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options);
  pet <- initVals$pet
  tdew <- initVals$tdew
  pva <- initVals$pva
  data <- initVals$data
  # /* estimate annual PET */
  ann_pet = pet * 365.25;
  
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
    # print(rmse_tdew)
    update_pva <- 1;
    tdew_save <- tdew
    initVals<-compute_srad_humidity_onetime_t1(ctrl, data, tdew, pva, ttmax0, flat_potrad, 
                                            slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options);
    pet <- initVals$pet
    tdew <- initVals$tdew
    pva <- initVals$pva
    data <- initVals$data
    rmse_tdew <- ((tdew-tdew_save)*(tdew-tdew_save))
    rmse_tdew <- rmse_tdew^0.5;
    iter <- iter + 1
  }
  
  # /* save humidity in output data structure */
  if (ctrl$outhum) {
    if (!ctrl$invp) {
      data$s_hum = pva
    }
  } else {
    # /* output humidity as vapor pressure deficit (Pa) */
    # /* calculate saturated VP at tday */
    # /* start vic_change */
    # /* pvs = 610.7 * exp(17.38 * data$s_tday/(239.0+data$s_tday)); */
    pvs = svp(data$s_tday);
    # /* end vic_change */
    vpd = pvs - pva;
    if (vpd < 0.0) vpd = 0.0;
    data$s_hum = vpd;
  }
  
  mt <- NULL
  mt <- data

  return(mt)
} # /* end of calc_srad_humidity_iterative() */

