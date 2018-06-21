# /* iterative estimation of shortwave radiation and humidity */
# /* Note: too many changes to maintain the start/end vic change comments */
# int calc_srad_humidity_iterative(const control_struct *ctrl,
#                                  const parameter_struct *p, data_struct *data,
#                                  double **tiny_radfract)
calc_tiny_radfract <- function(mt)
{
  tiny_radfract <- array(0, dim = c(366, 2880))
  daylength <- slope_potrad<- flat_potrad<- ttmax0 <- array(0, dim = c(366))
  
  p <- mt$p
  constants <- metGen$constants
  
  # /* optical airmass by degrees */
  optam <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
             4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
             8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
  #   /*****************************************
  #   *                                       *
  #   * start of the main radiation algorithm *
  #   *                                       *
  #   *****************************************/
  #   
  #   /* before starting the iterative algorithm between humidity and 
  # radiation, calculate all the variables that don't depend on 
  # humidity so they only get done once. */
  # 
  # /* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
  t1 = 1.0 - (constants$LR_STD * p$site_elev)/constants$T_STD;
  t2 = constants$G_STD / (constants$LR_STD * (constants$R/constants$MA));
  pratio = t1^t2
  
  # /* STEP (2) correct initial transmittance for elevation */ 
  trans1 = constants$TBASE^pratio
  
  # /* STEP (3) build 366-day array of ttmax0, potential rad, and daylength */
  
  # /* precalculate the transcendentals */
  lat = p$site_lat;
  # /* check for (+/-) 90 degrees latitude, throws off daylength calc */
  lat <- lat * constants$RADPERDEG;
  if (lat > 1.5707) lat = 1.5707;
  if (lat < -1.5707) lat = -1.5707;
  coslat = cos(lat);
  sinlat = sin(lat);
  cosslp = cos(p$site_slp * constants$RADPERDEG);
  sinslp = sin(p$site_slp * constants$RADPERDEG);
  cosasp = cos(p$site_asp * constants$RADPERDEG);
  sinasp = sin(p$site_asp * constants$RADPERDEG);
  # /* cosine of zenith angle for east and west horizons */
  coszeh = cos(1.570796 - (p$site_ehoriz * constants$RADPERDEG));
  coszwh = cos(1.570796 - (p$site_whoriz * constants$RADPERDEG));
  
  # /* sub-daily time and angular increment information */
  dt = constants$SRADDT;              #  /* set timestep */
  dh = dt / constants$SECPERRAD;      #  /* calculate hour-angle step */
  # /* start vic_change */
  tinystepspday = 86400/constants$SRADDT;
  # /* end vic_change */
  
  # /* begin loop through yeardays */
  for (i in 1:365) {
    # /* calculate cos and sin of declination */
    decl <- constants$MINDECL * cos(((i-1) + constants$DAYSOFF) * constants$RADPERDAY);
    cosdecl = cos(decl);
    sindecl = sin(decl);
    
    # /* do some precalculations for beam-slope geometry (bsg) */
    bsg1 = -sinslp * sinasp * cosdecl;
    bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl;
    bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl;
    
    # /* calculate daylength as a function of lat and decl */
    cosegeom = coslat * cosdecl;
    sinegeom = sinlat * sindecl;
    coshss = -(sinegeom) / cosegeom;
    if (coshss < -1.0) coshss = -1.0;   # /* 24-hr daylight */
    if (coshss > 1.0) coshss = 1.0;  #  /* 0-hr daylight */
    hss = acos(coshss);               # /* hour angle at sunset (radians) */
    # /* daylength (seconds) */
    daylength[i] = 2.0 * hss * constants$SECPERRAD;
    
    # /* start vic_change */
    if (daylength[i] > 86400) daylength[i] = 86400;
    # /* end vic_change */
    
    # /* solar constant as a function of yearday (W/m^2) */
    sc <- 1368.0 + 45.5*sin((2.0*pi*(i-1)/365.25) + 1.7);
    # /* extraterrestrial radiation perpendicular to beam, total over
    # the timestep (J) */
    dir_beam_topa = sc * dt;
    
    sum_trans = 0.0;
    sum_flat_potrad = 0.0;
    sum_slope_potrad = 0.0;
    
    # /* begin sub-daily hour-angle loop, from -hss to hss */
    for (h in seq(-hss,hss,dh)) {
      # /* precalculate cos and sin of hour angle */
      cosh = cos(h);
      sinh = sin(h);
      
      # /* calculate cosine of solar zenith angle */
      cza = cosegeom * cosh + sinegeom;
      
      # /* calculate cosine of beam-slope angle */
      cbsa = sinh * bsg1 + cosh * bsg2 + bsg3;
      
      # /* check if sun is above a flat horizon */
      if (cza > 0.0) {
        # /* when sun is above the ideal (flat) horizon, do all the
        # flat-surface calculations to determine daily total
        # transmittance, and save flat-surface potential radiation
        # for later calculations of diffuse radiation */
        
        # /* potential radiation for this time period, flat surface,
        # top of atmosphere */
        dir_flat_topa = dir_beam_topa * cza;
        
        # /* determine optical air mass */
        am = 1.0/(cza + 0.0000001);
        if (am > 2.9) {
          ami = (acos(cza)/constants$RADPERDEG) - 69;
          if (ami < 0)
            ami = 0;
          if (ami > 20)
            ami = 20;
          am = optam[ami+1];
        }
        
        # /* correct instantaneous transmittance for this optical
        # air mass */
        trans2 = trans1^am
        
        # /* instantaneous transmittance is weighted by potential
        # radiation for flat surface at top of atmosphere to get
        # daily total transmittance */
        sum_trans <- sum_trans + (trans2 * dir_flat_topa)
        
        # /* keep track of total potential radiation on a flat
        # surface for ideal horizons */
        sum_flat_potrad <- sum_flat_potrad + dir_flat_topa;
        
        # /* keep track of whether this time step contributes to
        # component 1 (direct on slope) */
        if ((h<0.0 && cza>coszeh && cbsa>0.0) ||  (h>=0.0 && cza>coszwh && cbsa>0.0)) {
          
          # /* sun between east and west horizons, and direct on
          # slope. this period contributes to component 1 */
          sum_slope_potrad <- sum_slope_potrad + ( dir_beam_topa * cbsa)
        }
        
      } else { # /* end if sun above ideal horizon */
        dir_flat_topa = -1; 
      }
      
      # /* start vic_change */
      tinystep = ((12L * 3600L + h * constants$SECPERRAD)/constants$SRADDT) + 1;
      if (tinystep < 1) tinystep <- 1;
      if (tinystep > (tinystepspday)) tinystep = tinystepspday;
      if (dir_flat_topa > 0) {
        tiny_radfract[i,tinystep] <- dir_flat_topa;
      } else {
        tiny_radfract[i,tinystep] <- 0;
      }
      # /* end vic_change */
      
    } # /* end of sub-daily hour-angle loop */
    
    # # /* start vic_change */
    if (daylength[i] && sum_flat_potrad > 0) {
      for (j in 1:tinystepspday) {
        tiny_radfract[i,j] <- tiny_radfract[i,j] / sum_flat_potrad;
      }
    }
    # # /* end vic_change */
    
    # # /* calculate maximum daily total transmittance and daylight average
    # # flux density for a flat surface and the slope */
    if (daylength[i]) {
      ttmax0[i] = sum_trans / sum_flat_potrad;
      # print(sum_trans)
      flat_potrad[i] = sum_flat_potrad / daylength[i];
      slope_potrad[i] = sum_slope_potrad / daylength[i];
    } else {
      ttmax0[i] = 0.0;
      flat_potrad[i] = 0.0;
      slope_potrad[i] = 0.0;
    }
    
  } # /* end of i=365 days loop */
  
  # /* force yearday 366 = yearday 365 */
  ttmax0[366] = ttmax0[365];
  flat_potrad[366] = flat_potrad[365];
  slope_potrad[366] = slope_potrad[365];
  daylength[366] = daylength[365];
  
  # /* start vic_change */
  for (j in 1:tinystepspday) tiny_radfract[366,j] = tiny_radfract[365,j];
  # /* end vic_change */
  
  mt$tiny_radfract <- tiny_radfract
  mt$ttmax0 <- ttmax0
  mt$flat_potrad <-flat_potrad
  mt$slope_potrad <- slope_potrad
  mt$daylength <- daylength
  
  return (mt) 
}  
