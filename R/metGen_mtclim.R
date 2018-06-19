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
  
  ndays <- ctrl$ndays;
  # /* calculate elevation difference in kilometers */
  dz <- (p$site_elev - p$base_elev)/1000.0;
  
  # /* apply lapse rate corrections to tmax and tmin */
  # /* Since tmax lapse rate usually has a larger absolute value than tmin
  # lapse rate, it is possible at high elevation sites for these corrections
  # to result in tmin > tmax. Check for that occurrence and force
  # tmin = corrected tmax - 0.5 deg C. */
  for (i in 1:ndays) {
    # /* lapse rate corrections */
    data$s_tmax[i] <- tmax <- data$tmax[i] + (dz * p$tmax_lr);
    data$s_tmin[i] <- tmin <- data$tmin[i] + (dz * p$tmin_lr);
    
    # /* derived temperatures */
    tmean = (tmax + tmin)/2.0;
    data$s_tday[i] = ((tmax - tmean)*constants$TDAYCOEF) + tmean;
  }
  
  mt$p <- p
  mt$ctrl <- ctrl
  mt$mtclim_data <- data
  
  return(mt)
}

# # /* calc_tair() calculates daily air temperatures */
# calc_tair_old<-function(tmin, tmax) {
#   # tmin <-4
#   # tmax<-5
#   # 
#   constants <- metGen$constants
#   
#   #   /* calculate elevation difference in kilometers */
#   #     dz = (p->site_elev - p->base_elev)/1000.0;
#   #   
#   #   /* apply lapse rate corrections to tmax and tmin */
#   #     /* Since tmax lapse rate usually has a larger absolute value than tmin
#   #   lapse rate, it is possible at high elevation sites for these corrections
#   #   to result in tmin > tmax. Check for that occurrence and force
#   #   tmin = corrected tmax - 0.5 deg C. */
#   #     for (i=0 ; i<ndays ; i++) {
#   #       /* lapse rate corrections */
#   #         data->s_tmax[i] = tmax = data->tmax[i] + (dz * p->tmax_lr);
#   #         data->s_tmin[i] = tmin = data->tmin[i] + (dz * p->tmin_lr);
#   
#   tmean <- (tmax + tmin)/2.0;
#   tair <- ((tmax - tmean)*constants$TDAYCOEF) + tmean;
#   
#   return(tair)  
# }

# /* calc_prcp() calculates daily total precipitation */
calc_prcp<-function(mt) {
  p <- mt$p
  ctrl <- mt$ctrl
  data <- mt$mtclim_data
  
  constants <- metGen$constants
  
  
  ndays <- ctrl$ndays;
  
  # /* start vic_change */
  ratio = -1.;
  if ( p$site_isoh < 1e-10 && p$base_isoh < 1e-10 ) {
    # /* If base_isoh and site_isoh are both small, set the ratio to 1.
    # This handles the case in which annual precip is 0, resulting in
    # base_isoh and site_isoh being 0 and their ratio being undefined. */
    ratio = 1.;
  }
  else if (p$base_isoh == 0) {
    stop("Error in calc_prcp(): base_isoh == 0 and site_isoh/base_isoh == NaN.");
  }
  else {
    ratio = p$site_isoh / p$base_isoh;
  }
  # /* end vic_change */
  
  for (i in 1:ndays) {
    data$s_prcp[i] = data$prcp[i] * ratio;
  }
  
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
  
  ndays <- ctrl$ndays;
  
  # /* first pass to initialize SWE array */
  snowpack = 0.0;
  for (i in 1:ndays)
  {
    newsnow = 0.0;
    snowmelt = 0.0;
    if (data$s_tmin[i] <= constants$SNOW_TCRIT) { newsnow = data$s_prcp[i];
    } else {
      snowmelt = constants$SNOW_TRATE * (data$s_tmin[i] - constants$SNOW_TCRIT);
    }
    snowpack <- snowpack + (newsnow - snowmelt)
    if (snowpack < 0.0) snowpack = 0.0;
    data$s_swe[i] = snowpack;
  }
  
  # /* use the first pass to set the initial snowpack conditions for the
  # first day of data */
  start_yday = data$yday[1];
  if (start_yday == 1) { prev_yday = 365 
  } else { prev_yday = start_yday-1; }
  count = 0;
  sum = 0.0;
  for (i in 2:ndays)
  {
    if (data$yday[i] == start_yday || data$yday[i] == prev_yday)
    {
      count <-count + 1;
      sum <- sum + data$s_swe[i];
    }
  }
  # /* Proceed with correction if there are valid days to reinitialize
  # the snowpack estiamtes. Otherwise use the first-pass estimate. */
  if (count)
  {
    snowpack = sum/count;
    for (i in 1:ndays)
    {
      newsnow = 0.0;
      snowmelt = 0.0;
      if (data$s_tmin[i] <= SNOW_TCRIT) {newsnow = data$s_prcp[i]; 
      } else { snowmelt = SNOW_TRATE * (data$s_tmin[i] - SNOW_TCRIT); }
      snowpack <-snowpack + (newsnow - snowmelt)
      if (snowpack < 0.0) snowpack = 0.0;
      data$s_swe[i] = snowpack;
    }
  }
  
  mt$p <- p
  mt$ctrl <- ctrl
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
      data$s_tskc[i] = (1.-data$s_tfmax[i]);
    } else {
      data$s_tskc[i] = sqrt((1.-data$s_tfmax[i])/0.65);
    }
    data$s_fdir[i] = pdir;
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
  
  
  
  return(initVals);
  
}

# /* iterative estimation of shortwave radiation and humidity */
# /* Note: too many changes to maintain the start/end vic change comments */
calc_srad_humidity_iterative_old <- function(elevation)
{
  tmin <- 4.6
  tmax <- 8.5
  prcp<- 1.4
  SW_PREC_THRESH <- 0
  ndays <- length(tmin)
  day <-100
  elevation<-434
  lat <- -8.25
  metGen$constants<-setConstants_new()
  constants <- metGen$constants
  tiny_radfract <- array(0, dim = c(2880))
  optam <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
             4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
             8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
  
  
  
  # /* calculate diurnal temperature range for transmittance calculations */
  if (tmax < tmin) 
  { dtr <- 0
  } else {
    dtr <- tmax-tmin;
  }
  
  # /* smooth dtr array: After Bristow and Campbell, 1984 */
  sm_dtr <- pulled_boxcar(dtr, 0)
  
  #     /*****************************************
  #     *                                       *
  #     * start of the main radiation algorithm *
  #     *                                       *
  #     *****************************************/
  #     
  #     /* before starting the iterative algorithm between humidity and 
  #   radiation, calculate all the variables that don't depend on 
  #   humidity so they only get done once. */
  #   
  #   /* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
  t1 <- 1.0 - (constants$LR_STD * elevation)/constants$T_STD;
  t2 <- constants$G_STD / (constants$LR_STD * (constants$R/constants$MA));
  
  pratio <-t1^t2
  
  #   /* STEP (2) correct initial transmittance for elevation */ 
  trans1 <- constants$TBASE^pratio
  
  #   /* STEP (3) build 366-day array of ttmax0, potential rad, and daylength */
  
  #   /* precalculate the transcendentals */
  #   /* check for (+/-) 90 degrees latitude, throws off daylength calc */
  lat <- lat * constants$RADPERDEG
  if (lat > 1.5707) lat <- 1.5707
  if (lat < -1.5707) lat <- -1.5707
  coslat <- cos(lat);
  sinlat <- sin(lat);
  asp<-slp<-0
  cosslp <- cos(slp * constants$RADPERDEG);
  sinslp <- sin(slp * constants$RADPERDEG);
  cosasp <- cos(asp * constants$RADPERDEG);
  sinasp <- sin(asp * constants$RADPERDEG);
  #   /* cosine of zenith angle for east and west horizons */
  ehoriz<- whoriz <- 0
  coszeh <- cos(1.570796 - (ehoriz * constants$RADPERDEG));
  coszwh <- cos(1.570796 - (whoriz * constants$RADPERDEG));
  
  #   /* sub-daily time and angular increment information */
  dt <- constants$SRADDT       #         /* set timestep */
  dh <- dt / constants$SECPERRAD;       # /* calculate hour-angle step */
  #   /* start vic_change */
  tinystepspday <- 86400/constants$SRADDT;
  #   /* end vic_change */
  
  #   /* begin loop through yeardays */
  # for (i in 1:365) {
  i<-day
  #   /* calculate cos and sin of declination */
  decl <- constants$MINDECL * cos(((i-1) + constants$DAYSOFF) * constants$RADPERDAY);
  cosdecl = cos(decl);
  sindecl = sin(decl);
  
  #   /* do some precalculations for beam-slope geometry (bsg) */
  bsg1 <- -sinslp * sinasp * cosdecl;
  bsg2 <- (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl;
  bsg3 <- (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl;
  
  #   /* calculate daylength as a function of lat and decl */
  cosegeom <- coslat * cosdecl;
  sinegeom <- sinlat * sindecl;
  coshss <- -(sinegeom) / cosegeom;
  if (coshss < -1.0)  coshss <- -1.0;  # /* 24-hr daylight */
  if (coshss > 1.0)   coshss <- 1.0;   # /* 0-hr daylight */
  hss <- acos(coshss);        #        /* hour angle at sunset (radians) */
  #   /* daylength (seconds) */
  daylength <- 2.0 * hss * constants$SECPERRAD
  
  #   /* start vic_change */
  if (daylength > 86400)  daylength[i] <- 86400;
  #   /* end vic_change */
  
  #   /* solar constant as a function of yearday (W/m^2) */
  sc <- 1368.0 + 45.5*sin((2.0*pi*(i-1)/365.25) + 1.7);
  #   /* extraterrestrial radiation perpendicular to beam, total over
  #   the timestep (J) */
  dir_beam_topa <- sc * dt;
  
  sum_trans <- 0.0;
  sum_flat_potrad <- 0.0;
  sum_slope_potrad <- 0.0;
  #   /* begin sub-daily hour-angle loop, from -hss to hss */
  for (h in seq(-hss,hss,dh)) {
    # h<-seq(-hss,hss,dh)[1]
    #   /* precalculate cos and sin of hour angle */
    cosh <- cos(h);
    sinh <- sin(h);
    
    #   /* calculate cosine of solar zenith angle */
    cza <- cosegeom * cosh + sinegeom;
    
    #   /* calculate cosine of beam-slope angle */
    cbsa <- sinh * bsg1 + cosh * bsg2 + bsg3;
    
    #   /* check if sun is above a flat horizon */
    if (cza > 0.0) {
      #   /* when sun is above the ideal (flat) horizon, do all the
      #   flat-surface calculations to determine daily total
      #   transmittance, and save flat-surface potential radiation
      #   for later calculations of diffuse radiation */
      #   
      #   /* potential radiation for this time period, flat surface,
      #   top of atmosphere */
      dir_flat_topa = dir_beam_topa * cza;
      
      #   /* determine optical air mass */
      am = 1.0/(cza + 0.0000001);
      if (am > 2.9) {
        ami = (acos(cza)/constants$RADPERDEG) - 69;
        if (ami < 0)
          ami = 0;
        if (ami > 20)
          ami = 20;
        am = optam[ami+1];
      }
      
      #   /* correct instantaneous transmittance for this optical
      # air mass */
      trans2 = trans1^am
      
      #   /* instantaneous transmittance is weighted by potential
      #   radiation for flat surface at top of atmosphere to get
      #   daily total transmittance */
      #   sum_trans += trans2 * dir_flat_topa;
      
      #   /* keep track of total potential radiation on a flat
      #   surface for ideal horizons */
      sum_flat_potrad <- sum_flat_potrad + dir_flat_topa;
      #   
      #   /* keep track of whether this time step contributes to
      #   component 1 (direct on slope) */
      if ((h<0.0 && cza>coszeh && cbsa>0.0) ||  (h>=0.0 && cza>coszwh && cbsa>0.0)) {
        
        #   /* sun between east and west horizons, and direct on
        #   slope. this period contributes to component 1 */
        sum_slope_potrad <- sum_slope_potrad + (dir_beam_topa * cbsa);
      }
    } #/* end if sun above ideal horizon */
    else { dir_flat_topa = -1;
    }
    #   /* start vic_change */
    tinystep = ((12L * 3600L + h * constants$SECPERRAD)/constants$SRADDT) + 1;
    if (tinystep < 1) tinystep <- 1;
    if (tinystep > (tinystepspday)) tinystep = tinystepspday;
    if (dir_flat_topa > 0) {
      tiny_radfract[tinystep] <- dir_flat_topa;
    } else {
      tiny_radfract[tinystep] <- 0;
    }
    #   /* end vic_change */
    
  } # /* end of sub-daily hour-angle loop */
  
  #   /* start vic_change */
  if (daylength && sum_flat_potrad > 0) {
    for (j in 1:tinystepspday) {
      tiny_radfract[j] <- tiny_radfract[j] / sum_flat_potrad
    }
  }
  
  #   /* end vic_change */
  
  #     /* calculate maximum daily total transmittance and daylight average
  #   flux density for a flat surface and the slope */
  if (daylength) {
    ttmax0 <- sum_trans / sum_flat_potrad;
    flat_potrad <- sum_flat_potrad / daylength;
    slope_potrad <- sum_slope_potrad / daylength;
  } else {
    ttmax0 <- 0.0;
    flat_potrad <- 0.0;
    slope_potrad <- 0.0;
  }
  #   
  # } /* end of i=365 days loop */
  
  # }
  # /* STEP (4)  calculate the sky proportion for diffuse radiation */
  #   /* uses the product of spherical cap defined by average horizon angle
  # and the great-circle truncation of a hemisphere. this factor does not
  # vary by yearday. */
  avg_horizon <- (ehoriz + whoriz)/2.0;
  horizon_scalar <- 1.0 - sin(avg_horizon * constants$RADPERDEG);
  if (slp > avg_horizon) {
    slope_excess = slp - avg_horizon;
  } else {
    slope_excess = 0.0;
  }
  
  if (2.0*avg_horizon > 180.0) {
    slope_scalar = 0.0;
  }else {
    slope_scalar = 1.0 - (slope_excess/(180.0 - 2.0*avg_horizon));
    if (slope_scalar < 0.0) 
      slope_scalar = 0.0;
  }
  sky_prop <- horizon_scalar * slope_scalar;
  
  # /* b parameter, and t_fmax not varying with Tdew, so these can be
  # calculated once, outside the iteration between radiation and humidity
  # estimates. Requires storing t_fmax in an array. */
  t_fmax <- array(NA, dim = ndays)
  s_tfmax <- array(NA, dim = ndays)
  for (i in 1:ndays) {	
    # /* b parameter from 30-day average of DTR */
    b <- constants$B0 + constants$B1 * exp(-constants$B2 * sm_dtr[i]);
    
    # /* proportion of daily maximum transmittance */
    t_fmax[i] = 1.0 - 0.9 * exp(-b * (dtr[i]^constants$C));
    
    # /* correct for precipitation if this is a rain day */
    if (prcp[i] > SW_PREC_THRESH) t_fmax[i] <- t_fmax[i]* constants$RAIN_SCALAR;
    s_tfmax[i] <- t_fmax[i];
  }
  
  
  # WF: move up
  metGen$internal$inVar$dewpt$supplied <- TRUE
  indewpt <- 10
  tdew <-NA
  
  # /* Initial values of vapor pressure, etc */
  if (metGen$internal$inVar$dewpt$supplied) {
    # /* Observed Tdew supplied */
    for (i in 1:ndays) {
      tdew[i] = indewpt[i];
    }
  } else {
    # /* Estimate Tdew */
    for (i in 1:ndays) {
      tdew[i] = s_tmin[i];
    }
  }
}

# #####################################################################################
# options <- NULL
# options$MTCLIM_SWE_CORR <- TRUE
# metGen$internal$inVar$insw$supplied <- TRUE    #ctrl->insw
# # data$s_srad
# # data$s_tfmax
# options$LW_CLOUD <- 1
# LW_CLOUD_DEARDORFF <- 1
# data<-NULL
# data$year <- NULL             # /* array of year values */
# data$yday <- NULL             # /* array of yearday values */
# data$tmax <- NULL          # /* array of base maximum temperature values */
# data$tmin <- NULL          # /* array of base minimum temperature values */
# data$prcp <- NULL          # /* array of base daily precipitation values */
# data$tdew <- NULL          # /* array of base dewpoint temperature values */
# data$s_tmax <- NULL        # /* array of site tmax values */
# data$s_tmin <- NULL        # /* array of site tmin values */
# data$s_tday <- NULL        # /* array of site daylight temperature values */
# data$s_prcp <- NULL        # /* array of site prcp values */
# data$s_hum <- NULL         # /* array of site humidity values (VPD or VP, Pa) */
# data$s_srad <- NULL        # /* array of site shortwave radiation values */
# data$s_dayl <- NULL        # /* array of site daylength values */
# data$s_swe <- NULL         # /* array of site snowpack values */
# # /* start vic_change */
# data$s_fdir <- NULL	 # /* array of site values of direct fraction of shortwave radiation */
# data$s_tskc <- NULL	 # /* array of cloudiness values */
# data$s_ppratio <- NULL # /* array of pet/prcp ratio values */
# data$s_ttmax <- NULL # /* array of clear sky transmittance values */
# data$s_tfmax <- NULL # /* array of cloud transmittance factor values */
# # /* end vic_change */
# 
# ndays <- metGen$derived$nday
# ctrl <- NULL ## NOT NEEDED?
# tdew <- 5
# pva <- 6
# ttmax0<-3
# flat_potrad <- NULL
#   slope_potrad <- NULL
# sky_prop<- NULL
# daylength<- NULL
# pet <- NULL
# parray<- NULL
# pa<- NULL
# dtr<- NULL
# 
# 
# 
# data$yday <- metGen$derived$inYDays

# # /* New function, not originally part of MTCLIM code */
# compute_srad_humidity_onetime <- function(ndays, ctrl, data, tdew, pva, 
#                                           ttmax0, flat_potrad, slope_potrad, 
#                                           sky_prop, daylength, pet, parray, 
#                                           pa, dtr) {
#   # (int ndays, const control_struct
#   # *ctrl, data_struct *data, double *tdew, double
#   # *pva, double *ttmax0, double *flat_potrad, double
#   # *slope_potrad, double sky_prop, double *daylength,
#   # double *pet, double *parray, double pa, double *dtr) {
# 
#   for (i in 1:ndays) {
#     yday <- data$yday[i]
# 
#     # /*** Compute SW radiation ***/
#     t_tmax <- ttmax0[yday] + ABASE * pva[i];
#     if (t_tmax < 0.0001) t_tmax <- 0.0001; # // this is mainly for the case of observed VP supplied, for which t_tmax sometimes ends up being negative (when potential radiation is low and VP is high)
#     data$s_ttmax[i] <- t_tmax;
# 
#     # /* final daily total transmittance */
#     t_final <- t_tmax * data$s_tfmax[i];
# 
#     # /* estimate fraction of radiation that is diffuse, on an
#     # instantaneous basis, from relationship with daily total
#     # transmittance in Jones (Plants and Microclimate, 1992)
#     # Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
#     # Fig 6.14, p. 122. */
#     pdif <- -1.25*t_final + 1.25;
#     if (pdif > 1.0) pdif <- 1.0;
#     if (pdif < 0.0) pdif <- 0.0;
# 
#     # /* estimate fraction of radiation that is direct, on an
#     # instantaneous basis */
#     pdir <- 1.0 - pdif;
# 
#     # /* the daily total radiation is estimated as the sum of the
#     # following two components:
#     #   1. The direct radiation arriving during the part of
#     # the day when there is direct beam on the slope.
#     # 2. The diffuse radiation arriving over the entire daylength
#     # (when sun is above ideal horizon).
#     # */
# 
#     # /* component 1 */
#     srad1 <- slope_potrad[yday] * t_final * pdir;
# 
#     # /* component 2 (diffuse) */
#     # /* includes the effect of surface albedo in raising the diffuse
#     # radiation for obstructed horizons */
#     srad2 <- flat_potrad[yday] * t_final * pdif * (sky_prop + constants$DIF_ALB*(1.0-sky_prop));
# 
#     # /* snow pack influence on radiation */
#     if (options$MTCLIM_SWE_CORR && data->s_swe[i] > 0.0) {
#       # /* snow correction in J/m2/day */
#       sc <- (1.32 + 0.096 * data$s_swe[i]) * 1e6;
#       # /* convert to W/m2 and check for zero daylength */
#       if (daylength[yday] > 0.0) { sc <- sc / daylength[yday];
#       } else { sc <- 0.0; }
#       # /* set a maximum correction of 100 W/m2 */
#       if (sc > 100.0) { sc <- 100.0; }
#     } else { sc <- 0.0; }
# 
#     # /* save daily radiation */
#     # /* save cloud transmittance when rad is an input */
#     if (metGen$internal$inVar$insw$supplied) {
#       potrad <- (srad1+srad2+sc)*daylength[yday]/t_final/86400;
#       if (potrad>0 && data$s_srad[i]>0 && daylength[yday]>0) {
#         data$s_tfmax[i] <- (data$s_srad[i])/(potrad*t_tmax); #//both of these are 24hr mean rad. here
#         if (data$s_tfmax[i] > 1.0) { data$s_tfmax[i] <- 1.0; }
#       } else {
#         data$s_tfmax[i] <- 1.0;
#       }
#     }
#     else {
#       data->s_srad[i] <- srad1 + srad2 +sc;
#     }
# 
#     # /* start vic_change */
#     if (options.LW_CLOUD == LW_CLOUD_DEARDORFF) {
#       data$s_tskc[i] <- (1.-data$s_tfmax[i]);
#     }
#     else {
#       data$s_tskc[i] <- sqrt((1.-data$s_tfmax[i])/0.65);
#     }
#     data$s_fdir[i] <- pdir;
#     # /* end vic_change */
# 
#   }
# 
#   # /*** Compute PET using SW radiation estimate, and update Tdew, pva ***/
#   for (i in 1:ndays) {
# 
#     tmink <- data$s_tmin[i] + KELVIN;
#     pet[i] <- calc_pet(data$s_srad[i],data$s_tday[i],pa,data$s_dayl[i]);
# 
#     # /* calculate ratio (PET/effann_prcp) and correct the dewpoint */
#     ratio <- pet[i]/parray[i];
#     data$s_ppratio[i] <- ratio*365.25;
#     ratio2 <- ratio*ratio;
#     ratio3 <- ratio2*ratio;
#     tdewk <- tmink*(-0.127 + 1.121*(1.003 - 1.444*ratio + 12.312*ratio2 - 32.766*ratio3) + 0.0006*(dtr[i]));
#     tdew[i] <- tdewk - KELVIN;
# 
#     # /* start vic_change */
#     # /* pva[i] = 610.7 * exp(17.38 * tdew[i] / (239.0 + tdew[i])); */
#     pva[i] <- svp(tdew[i]);
#     # /* end vic_change */
#   }
#   # return;
# }
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

