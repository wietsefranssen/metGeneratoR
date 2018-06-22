mtclim_init <- function(have_dewpt, have_shortwave, elevation, aspect, ehoriz, whoriz,
                        lat, prec, tmax, tmin, vp, yday, hourlyrad, 
                        tiny_radfract,
                        p, mtclim_data)
{
  # /* initialize the control structure */
  ctrl<-NULL
  # ctrl$ndays = metGen$derived$nday;
  
  ctrl$indewpt = 0;
  ctrl$invp = 0;
  if (have_dewpt) {
    if (have_dewpt == 1) {
      nrerror("have_dewpt not yet implemented for tdew; however you can supply observed vapor pressure and set have_dewpt to 2\n");
    }
    else if (have_dewpt == 2) {
      ctrl$invp = 1;
    }
  }
  if (have_shortwave) { ctrl$insw = 1;
  }
  else { ctrl$insw = 0;
  }
  ctrl$outhum = 1;	#	/* output vapor pressure */
  ctrl$inyear = 0;
  
  # /* initialize the parameter structure.  Meteorological variables are only
  #    calculated for the mean grid cell elevation.  The temperatures are lapsed
  #    outside of the mtclim code.  Therefore p$base_elev and p$site_elev are
  #    set to the same value.
  p <- NULL
  p$base_elev   = elevation;
  p$site_lat    = lat;
  p$site_elev   = elevation;
  p$site_slp    = 0 # slope;
  p$site_asp    = aspect;
  p$site_ehoriz = ehoriz;
  p$site_whoriz = whoriz;
  p$tmax_lr     = -1*metGen$constants$T_LAPSE*metGen$constants$METERS_PER_KM;	# /* not used since site_elev == base_elev */
  p$tmin_lr     = -1*metGen$constants$T_LAPSE*metGen$constants$METERS_PER_KM;	# /* not used since site_elev == base_elev */
  
  # # /* allocate space in the data arrays for input and output data */
  # if (data_alloc(ctrl, mtclim_data)) {
  #   stop("Error in data_alloc()... exiting\n");
  # }
  
  # /* initialize the data arrays with the vic input data */
  mtclim_data <- NULL
  mtclim_data$yday = yday
  mtclim_data$tmax = tmax;
  mtclim_data$tmin = tmin;
  if (ctrl$insw) {
    mtclim_data$s_srad = 0;
    for (j in 1:24) 
    {
      mtclim_data$s_srad <- mtclim_data$s_srad + hourlyrad[j];
    }
    mtclim_data$s_srad <- mtclim_data$s_srad / 24;
  }
  if (ctrl$invp) mtclim_data$s_hum = vp;
  # /* MTCLIM prcp in cm */
  mtclim_data$prcp = prec/10.; 
  if (have_dewpt==1)
    stop("have_dewpt not yet implemented ...\n");
  
  tinystepspday = 86400/metGen$constants$SRADDT;
  tiny_radfract <- array(0, dim = c(366,tinystepspday))
  
  mt<-NULL
  mt$mtclim_data<-mtclim_data
  mt$p <- p
  mt$ctrl <- ctrl
  mt$tiny_radfract<-tiny_radfract
  return(mt)
}
