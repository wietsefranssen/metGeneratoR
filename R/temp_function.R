# /* iterative estimation of shortwave radiation and humidity */
# /* Note: too many changes to maintain the start/end vic change comments */
# int calc_srad_humidity_iterative(const control_struct *ctrl,
#                                  const parameter_struct *p, data_struct *data,
#                                  double **tiny_radfract)

calc_rest <- function(mt, options)
{
  p <- mt$p
  ctrl <- mt$ctrl
  data <- mt$mtclim_data
  tiny_radfract <- mt$tiny_radfract
  ttmax0 <- mt$ttmax0
  flat_potrad <- mt$flat_potrad
  slope_potrad <- mt$slope_potrad
  daylength <- mt$daylength
  
  constants <- metGen$constants
  
  # /* optical airmass by degrees */
  optam <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
             4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
             8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
  
  # /* number of simulation days */
  ndays = ctrl$ndays;
  
  # /* local array memory allocation */
  dtr <- array(0, dim = c(ndays))
  parray <- array(0, dim = c(ndays))
  t_fmax <- array(0, dim = c(ndays))
  tdew <- array(0, dim = c(ndays))
  pet <- array(0, dim = c(ndays))
  pva <- array(0, dim = c(ndays))
  tdew_save <- array(0, dim = c(ndays))
  pva_save <- array(0, dim = c(ndays))
  window <- array(0, dim = c(ndays+90))
  
  # /* calculate diurnal temperature range for transmittance calculations */
  for (i in 1:ndays) {
    tmax = data$tmax[i];
    tmin = data$tmin[i];
    if (tmax < tmin) 
      tmax = tmin;
    dtr[i] = tmax-tmin;
  }
  
  # /* smooth dtr array: After Bristow and Campbell, 1984 */
  sm_dtr <- pulled_boxcar(dtr, 0)
  
  # /* calculate the annual total precip */
  sum_prcp = 0.0;
  for (i in 1:ndays) {
    sum_prcp<-sum_prcp + data$s_prcp[i];
  }
  ann_prcp = (sum_prcp/ndays) * 365.25;
  if (ann_prcp == 0.0) ann_prcp = 1.0;
  
  # /* Generate the effective annual precip, based on a 3-month
  # moving-window. Requires some special case handling for the
  # beginning of the record and for short records. */
  #   /* check if there are at least 90 days in this input file, if not,
  # use a simple total scaled to effective annual precip */
  if (ndays < 90) {
    sum_prcp = 0.0;
    for (i in 1:ndays) {
      sum_prcp <- sum_prcp + data$s_prcp[i];
    }
    effann_prcp = (sum_prcp/ndays) * 365.25;
    # /* if the effective annual precip for this period
    # is less than 8 cm, set the effective annual precip to 8 cm
    # to reflect an arid condition, while avoiding possible
    # division-by-zero errors and very large ratios (PET/Pann) */
    if (effann_prcp < 8.0) {
      effann_prcp = 8.0;
    }
    for (i in 1:ndays) {
      parray[i] = effann_prcp;
    }
  } else {
    # /* Check if the yeardays at beginning and the end of this input file
    # match up. If so, use parts of the three months at the end
    # of the input file to generate effective annual precip for
    # the first 3-months. Otherwise, duplicate the first 90 days
    # of the record. */
    start_yday = data$yday[1];
    end_yday = data$yday[ndays];
    if (start_yday != 1) {
      isloop = (end_yday == start_yday-1) ? 1 : 0;
    }
    else {
      isloop = (end_yday == 365 || end_yday == 366) ? 1 : 0;
    }
    
    # /* fill the first 90 days of window */
    for (i in 1:90) {
      if (isloop) {
        window[i] = data$s_prcp[ndays-90+i];  ## TODO!!!
      } else {
        window[i] = data$s_prcp[i];
      }
    }
    # /* fill the rest of the window array */
    for (i in 1:ndays) {
      window[i+90] = data$s_prcp[i];    ## TODO!!!
    }
    
    # /* for each day, calculate the effective annual precip from 
    # scaled 90-day total */
    for (i in 1:ndays) {
      sum_prcp = 0.0;
      for (j in 1:90) {
        sum_prcp <- sum_prcp + window[i+j];
      }
      sum_prcp = (sum_prcp/90.0) * 365.25;
      # /* if the effective annual precip for this 90-day period
      # is less than 8 cm, set the effective annual precip to 8 cm
      # to reflect an arid condition, while avoiding possible
      # division-by-zero errors and very large ratios (PET/Pann) */
      parray[i] = (sum_prcp < 8.0) ? 8.0 : sum_prcp;
    }
  } # /* end if ndays >= 90 */	
  
  # 
  # /* STEP (4)  calculate the sky proportion for diffuse radiation */
  # /* uses the product of spherical cap defined by average horizon angle
  # and the great-circle truncation of a hemisphere. this factor does not
  # vary by yearday. */
  avg_horizon = (p$site_ehoriz + p$site_whoriz)/2.0;
  horizon_scalar = 1.0 - sin(avg_horizon * constants$RADPERDEG);
  if (p$site_slp > avg_horizon) 
  {
    slope_excess = p$site_slp - avg_horizon;
  } else {
    slope_excess = 0.0;
  }
  if (2.0*avg_horizon > 180.0)
  {
    slope_scalar = 0.0;
  } else {
    slope_scalar = 1.0 - (slope_excess/(180.0 - 2.0*avg_horizon));
    if (slope_scalar < 0.0) slope_scalar = 0.0;
  }
  sky_prop = horizon_scalar * slope_scalar;
  
  # /* b parameter, and t_fmax not varying with Tdew, so these can be
  # calculated once, outside the iteration between radiation and humidity
  # estimates. Requires storing t_fmax in an array. */
  for (i in 1:ndays) {
    # /* b parameter from 30-day average of DTR */
    b = constants$B0 + constants$B1 * exp(-constants$B2 * sm_dtr[i]);
    
    # /* proportion of daily maximum transmittance */
    t_fmax[i] = 1.0 - 0.9 * exp(-b * dtr[i]^constants$C);
    
    # /* correct for precipitation if this is a rain day */
    if (data$prcp[i] > options$SW_PREC_THRESH) t_fmax[i] <- t_fmax[i] * constants$RAIN_SCALAR;
    data$s_tfmax[i] = t_fmax[i];
    
  }
  
  # /* Initial values of vapor pressure, etc */
  if (ctrl$indewpt) {
    # /* Observed Tdew supplied */
    for (i in 1:ndays) {
      tdew[i] = data$tdew[i];
    }
  } else {
    # /* Estimate Tdew */
    for (i in 1:ndays) {
      tdew[i] = data$s_tmin[i];
    }
  }
  if (ctrl$invp) {
    # /* Observed vapor pressure supplied */
    for (i in 1:ndays) {
      pva[i] = data$s_hum[i];
    }
  } else {
    # /* convert dewpoint to vapor pressure */
    for (i in 1:ndays) {
      # /* start vic_change */
      # /* pva[i] = 610.7 * exp(17.38 * tdew[i] / (239.0 + tdew[i])); */
      pva[i] = svp(tdew[i]);
      # /* end vic_change */
    }
  }
  
  # /* Other values needed for srad_humidity calculation */
  pa = atm_pres(p$site_elev);
  for (i in 1:ndays) {
    yday = data$yday[i];
    data$s_dayl[i] = daylength[yday];
    tdew_save[i] = tdew[i];
    pva_save[i] = pva[i];
  }
  
  # /* Initial estimates of solar radiation, cloud fraction, etc. */
  initVals <- compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, 
                                            slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options);
  pet <- initVals$pet
  tdew <- initVals$tdew
  pva <- initVals$pva
  data <- initVals$data
  # # /* estimate annual PET */
  sum_pet = 0.0;
  for (i in 1:ndays) {
    sum_pet <- sum_pet + pet[i];
  }
  ann_pet = (sum_pet/ndays) * 365.25;
  
  VP_ITER_ALWAYS <- 1
  VP_ITER_ANNUAL <- 2
  VP_ITER_CONVERGE <- 3
  
  # # /* Reset humidity terms if no iteration desired */
  if (ctrl$indewpt || ctrl$invp || (options$VP_ITER == VP_ITER_ANNUAL && ann_pet/ann_prcp >= 2.5) ) {
    for (i in 1:ndays) {
      tdew[i] = tdew_save[i];
      pva[i] = pva_save[i];
    }
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
    initVals<-compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options);
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
  mt$tiny_radfract <- tiny_radfract
  mt$ttmax0 <- ttmax0
  mt$flat_potrad <-flat_potrad
  mt$slope_potrad <- slope_potrad
  mt$daylength <- daylength
  
  return(mt)
} # /* end of calc_srad_humidity_iterative() */

calc_rest2 <- function(mt, options)
{
  p <- mt$p
  ctrl <- mt$ctrl
  data <- mt$mtclim_data
  tiny_radfract <- mt$tiny_radfract
  ttmax0 <- mt$ttmax0
  flat_potrad <- mt$flat_potrad
  slope_potrad <- mt$slope_potrad
  daylength <- mt$daylength
  
  constants <- metGen$constants
  
  # /* number of simulation days */
  ndays = ctrl$ndays;
  
  # /* local array memory allocation */
  dtr <- array(0, dim = c(ndays))
  parray <- array(0, dim = c(ndays))
  t_fmax <- array(0, dim = c(ndays))
  tdew <- array(0, dim = c(ndays))
  pet <- array(0, dim = c(ndays))
  pva <- array(0, dim = c(ndays))
  tdew_save <- array(0, dim = c(ndays))
  pva_save <- array(0, dim = c(ndays))
  window <- array(0, dim = c(ndays+90))
  
  # /* calculate diurnal temperature range for transmittance calculations */
  for (i in 1:ndays) {
    tmax = data$tmax[i];
    tmin = data$tmin[i];
    if (tmax < tmin) 
      tmax = tmin;
    dtr[i] = tmax-tmin;
  }
  
  # /* smooth dtr array: After Bristow and Campbell, 1984 */
  sm_dtr <- pulled_boxcar(dtr, 0)
  
  # /* calculate the annual total precip */
  sum_prcp = 0.0;
  for (i in 1:ndays) {
    sum_prcp<-sum_prcp + data$s_prcp[i];
  }
  ann_prcp = (sum_prcp/ndays) * 365.25;
  if (ann_prcp == 0.0) ann_prcp = 1.0;
  
  # /* Generate the effective annual precip, based on a 3-month
  # moving-window. Requires some special case handling for the
  # beginning of the record and for short records. */
  #   /* check if there are at least 90 days in this input file, if not,
  # use a simple total scaled to effective annual precip */
  if (ndays < 90) {
    sum_prcp = 0.0;
    for (i in 1:ndays) {
      sum_prcp <- sum_prcp + data$s_prcp[i];
    }
    effann_prcp = (sum_prcp/ndays) * 365.25;
    # /* if the effective annual precip for this period
    # is less than 8 cm, set the effective annual precip to 8 cm
    # to reflect an arid condition, while avoiding possible
    # division-by-zero errors and very large ratios (PET/Pann) */
    if (effann_prcp < 8.0) {
      effann_prcp = 8.0;
    }
    for (i in 1:ndays) {
      parray[i] = effann_prcp;
    }
  } else {
    # /* Check if the yeardays at beginning and the end of this input file
    # match up. If so, use parts of the three months at the end
    # of the input file to generate effective annual precip for
    # the first 3-months. Otherwise, duplicate the first 90 days
    # of the record. */
    start_yday = data$yday[0];
    end_yday = data$yday[ndays-1];
    if (start_yday != 1) {
      isloop = (end_yday == start_yday-1) ? 1 : 0;
    }
    else {
      isloop = (end_yday == 365 || end_yday == 366) ? 1 : 0;
    }
    
    # /* fill the first 90 days of window */
    for (i in 1:90) {
      if (isloop) {
        window[i] = data$s_prcp[ndays-90+i];  ## TODO!!!
      } else {
        window[i] = data$s_prcp[i];
      }
    }
    # /* fill the rest of the window array */
    for (i in 1:ndays) {
      window[i+90] = data$s_prcp[i];    ## TODO!!!
    }
    
    # /* for each day, calculate the effective annual precip from 
    # scaled 90-day total */
    for (i in 1:ndays) {
      sum_prcp = 0.0;
      for (j in 1:90) {
        sum_prcp <- sum_prcp + window[i+j];
      }
      sum_prcp = (sum_prcp/90.0) * 365.25;
      # /* if the effective annual precip for this 90-day period
      # is less than 8 cm, set the effective annual precip to 8 cm
      # to reflect an arid condition, while avoiding possible
      # division-by-zero errors and very large ratios (PET/Pann) */
      parray[i] = (sum_prcp < 8.0) ? 8.0 : sum_prcp;
    }
  } # /* end if ndays >= 90 */	
  
  # 
  # /* STEP (4)  calculate the sky proportion for diffuse radiation */
  # /* uses the product of spherical cap defined by average horizon angle
  # and the great-circle truncation of a hemisphere. this factor does not
  # vary by yearday. */
  avg_horizon = (p$site_ehoriz + p$site_whoriz)/2.0;
  horizon_scalar = 1.0 - sin(avg_horizon * constants$RADPERDEG);
  if (p$site_slp > avg_horizon) 
  {
    slope_excess = p$site_slp - avg_horizon;
  } else {
    slope_excess = 0.0;
  }
  if (2.0*avg_horizon > 180.0)
  {
    slope_scalar = 0.0;
  } else {
    slope_scalar = 1.0 - (slope_excess/(180.0 - 2.0*avg_horizon));
    if (slope_scalar < 0.0) slope_scalar = 0.0;
  }
  sky_prop <- horizon_scalar * slope_scalar;
  
  # /* b parameter, and t_fmax not varying with Tdew, so these can be
  # calculated once, outside the iteration between radiation and humidity
  # estimates. Requires storing t_fmax in an array. */
  for (i in 1:ndays) {
    # /* b parameter from 30-day average of DTR */
    b = constants$B0 + constants$B1 * exp(-constants$B2 * sm_dtr[i]);
    
    # /* proportion of daily maximum transmittance */
    t_fmax[i] = 1.0 - 0.9 * exp(-b * dtr[i]^constants$C);
    
    # /* correct for precipitation if this is a rain day */
    if (data$prcp[i] > options$SW_PREC_THRESH) t_fmax[i] <- t_fmax[i] * constants$RAIN_SCALAR;
    data$s_tfmax[i] = t_fmax[i];
    # print(t_fmax[i])
  }
  
  # /* Initial values of vapor pressure, etc */
  if (ctrl$indewpt) {
    # /* Observed Tdew supplied */
    for (i in 1:ndays) {
      tdew[i] = data$tdew[i];
    }
  } else {
    # /* Estimate Tdew */
    for (i in 1:ndays) {
      tdew[i] = data$s_tmin[i];
    }
  }
  if (ctrl$invp) {
    # /* Observed vapor pressure supplied */
    for (i in 1:ndays) {
      pva[i] = data$s_hum[i];
    }
  } else {
    # /* convert dewpoint to vapor pressure */
    for (i in 1:ndays) {
      # /* start vic_change */
      # /* pva[i] = 610.7 * exp(17.38 * tdew[i] / (239.0 + tdew[i])); */
      pva[i] = svp(tdew[i]);
      # /* end vic_change */
    }
  }
  
  # /* Other values needed for srad_humidity calculation */
  pa = atm_pres(p$site_elev);
  for (i in 1:ndays) {
    yday = data$yday[i];
    data$s_dayl[i] = daylength[yday];
    tdew_save[i] = tdew[i];
    pva_save[i] = pva[i];
  }
  
  # /* Initial estimates of solar radiation, cloud fraction, etc. */
  initVals <- compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, 
                                            slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options);
  pet <- initVals$pet
  tdew <- initVals$tdew
  pva <- initVals$pva
  data <- initVals$data
  # # /* estimate annual PET */
  sum_pet = 0.0;
  for (i in 1:ndays) {
    sum_pet <- sum_pet + pet[i];
  }
  ann_pet = (sum_pet/ndays) * 365.25;
  
  VP_ITER_ALWAYS <- 1
  VP_ITER_ANNUAL <- 2
  VP_ITER_CONVERGE <- 3
  
  # # /* Reset humidity terms if no iteration desired */
  if (ctrl$indewpt || ctrl$invp || (options$VP_ITER == VP_ITER_ANNUAL && ann_pet/ann_prcp >= 2.5) ) {
    for (i in 1:ndays) {
      tdew[i] = tdew_save[i];
      pva[i] = pva_save[i];
    }
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
    initVals<-compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options);
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
  mt$tiny_radfract <- tiny_radfract
  mt$ttmax0 <- ttmax0
  mt$flat_potrad <-flat_potrad
  mt$slope_potrad <- slope_potrad
  mt$daylength <- daylength
  
  return(mt)
} # /* end of calc_srad_humidity_iterative() */

