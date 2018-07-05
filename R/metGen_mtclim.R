# /******************************************************************************
#   * @brief        This routine computes the saturated vapor pressure
# *
#   * @note         Handbook of Hydrology eqn 4.2.2.
# ******************************************************************************/
svp <-function(temp)
{
  # // Saturation Vapor Pressure Parameters
  param_SVP_A = 0.61078;
  param_SVP_B = 17.269;
  param_SVP_C = 237.3;
  SVP = param_SVP_A * exp((param_SVP_B * temp) / (param_SVP_C + temp));
  
  if (temp < 0) {
    SVP <-SVP  *( 1.0 + .00972 * temp + .000042 * temp * temp)
  }
  
  return (SVP * metGen$constants$PA_PER_KPA);
}

# /* atm_pres() calculates the atmospheric pressure as a function of elevation */
atm_pres <- function(elev)
{
  # /* daily atmospheric pressure (Pa) as a function of elevation (m) */
  # /* From the discussion on atmospheric statics in:
  # Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics, 2nd
  # Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
  # (p. 168)
  # */
  
  t1 = 1.0 - (metGen$constants$LR_STD * elev)/metGen$constants$T_STD;
  t2 = metGen$constants$G_STD / (metGen$constants$LR_STD * (metGen$constants$R / metGen$constants$MA));
  pa = metGen$constants$P_STD * t1^t2
  
  return(pa);
}

# /* calc_tair() calculates daily air temperatures */
calc_tair<-function(mt) {
  p <- mt$p
  ctrl <- mt$ctrl
  data <- mt$mtclim_data
  
  constants <- metGen$constants
  
  # /* calculate elevation difference in kilometers */
  dz <- (p$site_elev - p$base_elev)/1000.0;
  
  # /* apply lapse rate corrections to tmax and tmin */
  # /* Since tmax lapse rate usually has a larger absolute value than tmin
  # lapse rate, it is possible at high elevation sites for these corrections
  # to result in tmin > tmax. Check for that occurrence and force
  # tmin = corrected tmax - 0.5 deg C. */
  # /* lapse rate corrections */
  data$s_tmax <- tmax <- data$tmax + (dz * p$tmax_lr);
  data$s_tmin <- tmin <- data$tmin + (dz * p$tmin_lr);
  
  # /* derived temperatures */
  tmean = (tmax + tmin)/2.0;
  data$s_tday = ((tmax - tmean)*constants$TDAYCOEF) + tmean;
  
  mt$p <- p
  mt$ctrl <- ctrl
  mt$mtclim_data <- data
  
  return(mt)
}

# /* calc_prcp() calculates daily total precipitation */
calc_prcp<-function(mt) {
  p <- mt$p
  ctrl <- mt$ctrl
  data <- mt$mtclim_data
  
  constants <- metGen$constants
  
  
  ndays <- ctrl$ndays;
  
  data$s_prcp = data$prcp
  
  mt$p <- p
  mt$ctrl <- ctrl
  mt$mtclim_data <- data
  
  return(mt)
}

# /* snowpack() estimates the accumulation and melt of snow for radiation
# algorithm corrections */
snowpack<-function(mt) {
  p <- mt$p
  ctrl <- mt$ctrl
  data <- mt$mtclim_data
  constants <- metGen$constants
  
  # /* first pass to initialize SWE array */
  snowpack = 0.0;
  newsnow = 0.0;
  snowmelt = 0.0;
  if (data$s_tmin <= constants$SNOW_TCRIT) { newsnow = data$s_prcp;
  } else {
    snowmelt = constants$SNOW_TRATE * (data$s_tmin - constants$SNOW_TCRIT);
  }
  snowpack <- snowpack + (newsnow - snowmelt)
  if (snowpack < 0.0) snowpack = 0.0;
  data$s_swe = snowpack;
  
  
  # /* use the first pass to set the initial snowpack conditions for the
  # first day of data */
  start_yday = data$yday[1];
  if (start_yday == 1) { prev_yday = 365 
  } else { prev_yday = start_yday-1; }
  count = 0;
  sum = 0.0;
  
  # /* Proceed with correction if there are valid days to reinitialize
  # the snowpack estiamtes. Otherwise use the first-pass estimate. */
  if (count)
  {
    snowpack = sum/count;
    newsnow = 0.0;
    snowmelt = 0.0;
    if (data$s_tmin <= constants$SNOW_TCRIT) {newsnow = data$s_prcp; 
    } else { snowmelt = constants$SNOW_TRATE * (data$s_tmin - constants$SNOW_TCRIT); }
    snowpack <-snowpack + (newsnow - snowmelt)
    if (snowpack < 0.0) snowpack = 0.0;
    data$s_swe = snowpack;
  }
  
  mt$p <- p
  mt$mtclim_data <- data
  
  return(mt)
}

# /* calc_pet() calculates the potential evapotranspiration for aridity 
# corrections in calc_vpd(), according to Kimball et al., 1997 */
calc_pet <- function(rad, ta, pa, dayl)
{
  # /* input parameters and units :
  #   double rad      (W/m2)  daylight average incident shortwave radiation
  # double ta       (deg C) daylight average air temperature
  # double pa       (Pa)    air pressure
  # double dayl     (s)     daylength
  # */
  #
  #   double rnet;       /* (W m-2) absorbed shortwave radiation avail. for ET */
  #   double lhvap;      /* (J kg-1) latent heat of vaporization of water */
  #   double gamma;      /* (Pa K-1) psychrometer parameter */
  #   double dt = 0.2;   /* offset for saturation vapor pressure calculation */
  #   double t1, t2;     /* (deg C) air temperatures */
  #   double pvs1, pvs2; /* (Pa)   saturated vapor pressures */
  #   double pet;        /* (kg m-2 day-1) potential evapotranspiration */
  #   double s;          /* (Pa K-1) slope of saturated vapor pressure curve */
  dt <- 0.3
  # /* calculate absorbed radiation, assuming albedo = 0.2  and ground
  # heat flux = 10% of absorbed radiation during daylight */
  rnet = rad * 0.72;
  
  # /* calculate latent heat of vaporization as a function of ta */
  lhvap = 2.5023e6 - 2430.54 * ta;
  
  # /* calculate the psychrometer parameter: gamma = (cp pa)/(lhvap epsilon)
  # where:
  #   cp       (J/kg K)   specific heat of air
  # epsilon  (unitless) ratio of molecular weights of water and air
  # */
  gamma = metGen$constants$CP * pa / (lhvap * metGen$constants$EPS);
  
  # /* estimate the slope of the saturation vapor pressure curve at ta */
  # /* temperature offsets for slope estimate */
  t1 = ta+dt;
  t2 = ta-dt;
  
  # /* calculate saturation vapor pressures at t1 and t2, using formula from
  # Abbott, P.F., and R.C. Tabony, 1985. The estimation of humidity parameters.
  # Meteorol. Mag., 114:49-56.
  # */
  #   /* start vic_change */
  #   /* pvs1 = 610.7 * exp(17.38 * t1 / (239.0 + t1)); */
  pvs1 = svp(t1);
  # /* pvs2 = 610.7 * exp(17.38 * t2 / (239.0 + t2)); */
  pvs2 = svp(t2);
  # /* end vic_change */
  
  # /* calculate slope of pvs vs. T curve near ta */
  s = (pvs1-pvs2) / (t1-t2);
  
  # /* calculate PET using Priestly-Taylor approximation, with coefficient
  # set at 1.26. Units of result are kg/m^2/day, equivalent to mm water/day */
  pet = (1.26 * (s/(s+gamma)) * rnet * dayl)/lhvap;
  
  # /* return a value in centimeters/day, because this value is used in a ratio
  # to annual total precip, and precip units are centimeters */
  return (pet/10.0);
}

# # /* calc_prcp() calculates daily total precipitation */
# calc_prcp_old<-function(prcp) {
#   
#   return(prcp)  
# }
# /* New function, not originally part of MTCLIM code */
# void compute_srad_humidity_onetime(int ndays, const control_struct *ctrl, data_struct *data, double *tdew, double *pva, double *ttmax0, double *flat_potrad, double *slope_potrad, double sky_prop, double *daylength, double *pet, double *parray, double pa, double *dtr) {
compute_srad_humidity_onetime <- function(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, 
                                          slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options) {
  for (i in 1:ndays) {
    yday = data$yday[i];
    
    #   /*** Compute SW radiation ***/
    t_tmax = ttmax0[yday] + metGen$constants$ABASE * pva[i];
    if (t_tmax < 0.0001) t_tmax = 0.0001; # // this is mainly for the case of observed VP supplied, for which t_tmax sometimes ends up being negative (when potential radiation is low and VP is high)
    data$s_ttmax[i] = t_tmax;
    
    #     /* final daily total transmittance */
    t_final = t_tmax * data$s_tfmax[i];
    #     /* estimate fraction of radiation that is diffuse, on an
    #     instantaneous basis, from relationship with daily total
    #     transmittance in Jones (Plants and Microclimate, 1992)
    #     Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
    #     Fig 6.14, p. 122. */
    pdif = -1.25*t_final + 1.25;
    if (pdif > 1.0) pdif = 1.0;
    if (pdif < 0.0) pdif = 0.0;
    
    #     /* estimate fraction of radiation that is direct, on an
    #     instantaneous basis */
    pdir = 1.0 - pdif;
    
    #     /* the daily total radiation is estimated as the sum of the
    #     following two components:
    #       1. The direct radiation arriving during the part of
    #     the day when there is direct beam on the slope.
    #     2. The diffuse radiation arriving over the entire daylength
    #     (when sun is above ideal horizon).
    #     */
    #       
    #       /* component 1 */
    srad1 = slope_potrad[yday] * t_final * pdir;
    
    #     /* component 2 (diffuse) */
    #       /* includes the effect of surface albedo in raising the diffuse
    #     radiation for obstructed horizons */
    srad2 = flat_potrad[yday] * t_final * pdif * (sky_prop + metGen$constants$DIF_ALB*(1.0-sky_prop));
    
    #     /* snow pack influence on radiation */
    if (options$MTCLIM_SWE_CORR && data$s_swe[i] > 0.0) {
      #         /* snow correction in J/m2/day */
      sc = (1.32 + 0.096 * data$s_swe[i]) * 1e6;
      #           /* convert to W/m2 and check for zero daylength */
      if (daylength[yday] >  0.0) {
        sc<-sc / daylength[yday] 
      } else { 
        sc = 0.0 
      }
      #             /* set a maximum correction of 100 W/m2 */
      if (sc > 100.0) sc = 100.0;
    } else { sc = 0.0 }
    
    #     /* save daily radiation */
    #       /* save cloud transmittance when rad is an input */
    if (ctrl$insw) {
      
      potrad = (srad1+srad2+sc)*daylength[yday]/t_final/86400;
      if (potrad>0 && data$s_srad[i]>0 && daylength[yday]>0) {
        data$s_tfmax[i] = (data$s_srad[i])/(potrad*t_tmax); #//both of these are 24hr mean rad. here
        # print(data$s_tfmax[i])
        if (data$s_tfmax[i] > 1.0) data$s_tfmax[i] = 1.0 
      } else {
        data$s_tfmax[i] = 1.0;
      }
    } else {
      data$s_srad[i] = srad1 + srad2 +sc;
    }
    
    #     /* start vic_change */
    LW_CLOUD_DEARDORFF <- 1
    if (options$LW_CLOUD == LW_CLOUD_DEARDORFF) {
      data$s_tskc[i] = (1. - data$s_tfmax[i]);
    } else {
      data$s_tskc[i] = sqrt((1. - data$s_tfmax[i])/0.65);
    }
    data$s_fdir[i] = pdir;
    # print(pdir)
    #     /* end vic_change */
    
  }
  
  KELVIN<- 273.15
  # /*** Compute PET using SW radiation estimate, and update Tdew, pva ***/
  for (i in 1:ndays) {
    tmink = data$s_tmin[i] + KELVIN;
    pet[i] = calc_pet(data$s_srad[i],data$s_tday[i],pa,data$s_dayl[i]);
    
    #     /* calculate ratio (PET/effann_prcp) and correct the dewpoint */
    ratio = pet[i]/parray[i];
    data$s_ppratio[i] = ratio*365.25;
    ratio2 = ratio*ratio;
    ratio3 = ratio2*ratio;
    tdewk = tmink*(-0.127 + 1.121*(1.003 - 1.444*ratio + 12.312*ratio2
                                   - 32.766*ratio3) + 0.0006*(dtr[i]));
    tdew[i] = tdewk - KELVIN;
    
    #       /* start vic_change */
    #         /* pva[i] = 610.7 * exp(17.38 * tdew[i] / (239.0 + tdew[i])); */
    pva[i] = svp(tdew[i]);
    # /* end vic_change */
    
  }
  
  initVals <- NULL
  initVals$pet <- pet
  initVals$pva <- pva
  initVals$tdew <- tdew
  initVals$data <- data
  
  
  # print(initVals)
  
  return(initVals);
  
}


#####################################################################################

# /* pulled_boxcar() calculates a moving average of antecedent values in an
# array, using either a ramped (w_flag=1) or a flat (w_flag=0) weighting */	
# int pulled_boxcar(double *input,double *output,int n,int w,int w_flag)
pulled_boxcar <- function(input,  w_flag) {
  # input <- c(1:35)
  # w_flag <- 0
  
  n <- length(input)
  w <- 30
  
  if (w > n) {
    w <- n
  }
  
  wt <- array(NA, dim = c(w))
  output <- array(0, dim = c(n))
  
  # /* when w_flag != 0, use linear ramp to weight tails,
  # otherwise use constant weight */
  sum_wt <- 0.0;
  if (w_flag) {
    for (i in 1:w) {
      wt[i] = i+1
      sum_wt <- sum_wt + wt[i]
    }
  } else {
    for (i in 1:w) { 	
      wt[i] <- 1.0
      sum_wt <- sum_wt + wt[i]
    }
  }
  
  # /* fill the output array, starting with the point where a full
  # boxcar can be calculated */
  for (i in w:n) {
    # i<-4
    total <- 0.0
    for (j in 1:w) {
      total <- total + input[i-(w+1)+j+1] * wt[j];
    }
    output[i] = total/sum_wt;
  }
  
  # /* fill the first w elements of the output array with the value from
  # the first full boxcar */
  for (i in 1:(w-1)) {
    output[i] = output[w];
  }
  return(output)
}
# /* end of pulled_boxcar() */  
# input <- c(10.96,10.96,12.05, 7.30)
# # input <- c(12.08,12.08,12.60, 8.30)
# output <- pulled_boxcar(input,  0)
# print(output)

