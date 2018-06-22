# /* New function, not originally part of MTCLIM code */
# void compute_srad_humidity_onetime(int ndays, const control_struct *ctrl, data_struct *data, double *tdew, double *pva, double *ttmax0, double *flat_potrad, double *slope_potrad, double sky_prop, double *daylength, double *pet, double *parray, double pa, double *dtr) {
compute_srad_humidity_onetime_t1 <- function(ctrl, data, tdew, pva, ttmax0, flat_potrad, 
                                          slope_potrad, sky_prop, daylength, pet, parray, pa, dtr, options) {
    yday = data$yday
    sky_prop <- 1
    #   /*** Compute SW radiation ***/
    t_tmax = ttmax0 + metGen$constants$ABASE * pva;
    if (t_tmax < 0.0001) t_tmax = 0.0001; # // this is mainly for the case of observed VP supplied, for which t_tmax sometimes ends up being negative (when potential radiation is low and VP is high)
    data$s_ttmax = t_tmax;
    
    #     /* final daily total transmittance */
    t_final = t_tmax * data$s_tfmax;
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
    srad1 = slope_potrad * t_final * pdir;
    
    #     /* component 2 (diffuse) */
    #       /* includes the effect of surface albedo in raising the diffuse
    #     radiation for obstructed horizons */
    srad2 = flat_potrad * t_final * pdif * (sky_prop + metGen$constants$DIF_ALB*(1.0-sky_prop));
    
    #     /* snow pack influence on radiation */
    if (options$MTCLIM_SWE_CORR && data$s_swe > 0.0) {
      #         /* snow correction in J/m2/day */
      sc = (1.32 + 0.096 * data$s_swe) * 1e6;
      #           /* convert to W/m2 and check for zero daylength */
      if (daylength >  0.0) {
        sc<-sc / daylength 
      } else { 
        sc = 0.0 
      }
      #             /* set a maximum correction of 100 W/m2 */
      if (sc > 100.0) sc = 100.0;
    } else { sc = 0.0 }
    
    #     /* save daily radiation */
    #       /* save cloud transmittance when rad is an input */
    if (ctrl$insw) {
      
      potrad = (srad1+srad2+sc)*daylength/t_final/86400;
      if (potrad>0 && data$s_srad>0 && daylength>0) {
        data$s_tfmax = (data$s_srad)/(potrad*t_tmax); #//both of these are 24hr mean rad. here
        # print(data$s_tfmax)
        if (data$s_tfmax > 1.0) data$s_tfmax = 1.0 
      } else {
        data$s_tfmax = 1.0;
      }
    } else {
      data$s_srad = srad1 + srad2 +sc;
    }
    
    #     /* start vic_change */
    LW_CLOUD_DEARDORFF <- 1
    if (options$LW_CLOUD == LW_CLOUD_DEARDORFF) {
      data$s_tskc = (1. - data$s_tfmax);
    } else {
      data$s_tskc = sqrt((1. - data$s_tfmax)/0.65);
    }
    data$s_fdir = pdir;
    # print(pdir)
    #     /* end vic_change */
    
  
  
  KELVIN<- 273.15
  # /*** Compute PET using SW radiation estimate, and update Tdew, pva ***/
    tmink = data$s_tmin + KELVIN;
    pet = calc_pet(data$s_srad,data$s_tday,pa,data$s_dayl);
    
    #     /* calculate ratio (PET/effann_prcp) and correct the dewpoint */
    ratio = pet/parray;
    data$s_ppratio = ratio*365.25;
    ratio2 = ratio*ratio;
    ratio3 = ratio2*ratio;
    tdewk = tmink*(-0.127 + 1.121*(1.003 - 1.444*ratio + 12.312*ratio2
                                   - 32.766*ratio3) + 0.0006*(dtr));
    tdew = tdewk - KELVIN;
    
    #       /* start vic_change */
    #         /* pva = 610.7 * exp(17.38 * tdew / (239.0 + tdew)); */
    pva = svp(tdew);
    # /* end vic_change */
  
  initVals <- NULL
  initVals$pet <- pet
  initVals$pva <- pva
  initVals$tdew <- tdew
  initVals$data <- data
  
  return(initVals);
  
}
