calc_tiny_radfract_rest <- function(mt, options)
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
  
  #############################################################################################################
  data <- mt$mtclim_data
  
  # /* calculate diurnal temperature range for transmittance calculations */
  tmax = data$tmax;
  tmin = data$tmin;
  if (tmax < tmin) 
    tmax = tmin;
  dtr = tmax-tmin;
  
  # /* smooth dtr array: After Bristow and Campbell, 1984 */
  sm_dtr <- pulled_boxcar(dtr, 0)
  
  # /* calculate the annual total precip */
  sum_prcp = 0.0;
  sum_prcp<-sum_prcp + data$s_prcp;
  ann_prcp = sum_prcp * 365.25;
  if (ann_prcp == 0.0) ann_prcp = 1.0;
  
  # /* Generate the effective annual precip, based on a 3-month
  # moving-window. Requires some special case handling for the
  # beginning of the record and for short records. */
  #   /* check if there are at least 90 days in this input file, if not,
  # use a simple total scaled to effective annual precip */
  sum_prcp = 0.0;
  sum_prcp <- sum_prcp + data$s_prcp;
  effann_prcp = sum_prcp * 365.25;
  # /* if the effective annual precip for this period
  # is less than 8 cm, set the effective annual precip to 8 cm
  # to reflect an arid condition, while avoiding possible
  # division-by-zero errors and very large ratios (PET/Pann) */
  if (effann_prcp < 8.0) {
    effann_prcp = 8.0;
  }
  parray = effann_prcp;
  
  #   /*****************************************
  #   *                                       *
  #   *                                       *
  #   *****************************************/
  
  
  #####################
  # /* calculate diurnal temperature range for transmittance calculations */
  tmax = data$tmax;
  tmin = data$tmin;
  if (tmax < tmin) 
    tmax = tmin;
  dtr = tmax-tmin;
  
  # /* smooth dtr array: After Bristow and Campbell, 1984 */
  sm_dtr <- pulled_boxcar(dtr, 0)
  
  # /* calculate the annual total precip */
  sum_prcp = 0.0;
  sum_prcp<-sum_prcp + data$s_prcp;
  ann_prcp = sum_prcp * 365.25;
  if (ann_prcp == 0.0) ann_prcp = 1.0;
  
  # /* Generate the effective annual precip, based on a 3-month
  # moving-window. Requires some special case handling for the
  # beginning of the record and for short records. */
  #   /* check if there are at least 90 days in this input file, if not,
  # use a simple total scaled to effective annual precip */
  sum_prcp = 0.0;
  sum_prcp <- sum_prcp + data$s_prcp;
  effann_prcp = sum_prcp * 365.25;
  # /* if the effective annual precip for this period
  # is less than 8 cm, set the effective annual precip to 8 cm
  # to reflect an arid condition, while avoiding possible
  # division-by-zero errors and very large ratios (PET/Pann) */
  if (effann_prcp < 8.0) {
    effann_prcp = 8.0;
  }
  parray = effann_prcp;
  
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
  # /* b parameter from 30-day average of DTR */
  b = constants$B0 + constants$B1 * exp(-constants$B2 * sm_dtr);
  
  # /* proportion of daily maximum transmittance */
  t_fmax = 1.0 - 0.9 * exp(-b * dtr^constants$C);
  
  # /* correct for precipitation if this is a rain day */
  if (data$prcp > options$SW_PREC_THRESH) t_fmax <- t_fmax * constants$RAIN_SCALAR;
  data$s_tfmax = t_fmax;
  
  #####################
  mt$p <- p
  mt$mtclim_data <- data
  mt$tiny_radfract <- tiny_radfract
  mt$ttmax0 <- ttmax0
  mt$flat_potrad <-flat_potrad
  mt$slope_potrad <- slope_potrad
  mt$daylength <- daylength
  mt$sky_prop <- sky_prop 
  return(mt)
} # /* end of calc_srad_humidity_iterative() */

