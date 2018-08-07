setConstants <- function() {
  constants <- NULL
  # /* parameters for the Tair algorithm */
  constants$TDAYCOEF    <- 0.45 #  /* (dim) daylight air temperature coefficient (dim) */
  
  # # /* parameters for the snowpack algorithm */
  constants$SNOW_TCRIT  <- -6.0 #  /* (deg C) critical temperature for snowmelt */
  constants$SNOW_TRATE  <- 0.042 # /* (cm/degC/day) snowmelt rate */
  # 
  # # /* parameters for the radiation algorithm */
  constants$TBASE   <-    0.870 # /* (dim) max inst. trans., 0m, nadir, dry atm */
  constants$ABASE   <-  -6.1e-5 # /* (1/Pa) vapor pressure effect on transmittance */
  constants$C       <-      1.5 # /* (dim) radiation parameter */
  constants$B0      <-    0.031 # /* (dim) radiation parameter */
  constants$B1      <-    0.201 # /* (dim) radiation parameter */
  constants$B2      <-    0.185 # /* (dim) radiation parameter */
  constants$RAIN_SCALAR <- 0.75 # /* (dim) correction to trans. for rain day */
  constants$DIF_ALB     <-  0.6 # /* (dim) diffuse albedo for horizon correction */
  # constants$SC_INT      <- 1.32 # /* (MJ/m2/day) snow correction intercept */
  # constants$SC_SLOPE  <-  0.096 # /* (MJ/m2/day/cm) snow correction slope */
  # 
  constants$TDAYCOEF <- 0.45  #   /* daylight air temperature coefficient (dim) */
  
  constants$SECPERRAD <- 13750.9871   #  /* seconds per radian of hour angle */
  constants$RADPERDAY <- 0.017214    #   /* radians of Earth orbit per julian day */
  constants$RADPERDEG  <- 0.01745329  #   /* radians per degree */
  constants$MINDECL <-  -0.4092797   #    /* minimum declination (radians) */
  constants$DAYSOFF <-  11.25       #     /* julian day offset of winter solstice */
  # /* start vic_change */
  constants$SRADDT <- 30.0       #      /* timestep for radiation routine (seconds) */
  #   /* Note:  Make sure that 3600 % SRADDT == 0 */
  #   /* end vic_change */
  #   
  constants$MA   <-    28.9644e-3    #  /* (kg mol-1) molecular weight of air */
  #   constants$MW       18.0148e-3  #    /* (kg mol-1) molecular weight of water */
  constants$R      <-  8.3143       #   /* (m3 Pa mol-1 K-1) gas law constant */
  constants$G_STD  <-  9.80665       #  /* (m s-2) standard gravitational accel. */ 
  constants$P_STD  <-  101325.0   #     /* (Pa) standard pressure at 0.0 m elevation */
  constants$T_STD   <- 288.15      #    /* (K) standard temp at 0.0 m elevation  */
    constants$CP   <-    1010.0    #      /* (J kg-1 K-1) specific heat of air */
  constants$LR_STD <-  0.0065      #    /* (-K m-1) standard temperature lapse rate */
  constants$T_LAPSE   <-   -0.0065   # /* temperature lapse rate of US Std Atmos in C/km */
  constants$METERS_PER_KM <- 1000. #
  constants$PA_PER_KPA <- 1000 # /**< Pa per kPa */
  constants$EPS <- 0.62196351
  
  ## New  
  constants$DEG_PER_REV <- 360.0       # Number of degrees in full revolution
  constants$SEC_PER_RAD <- 13750.9871  # seconds per radian of hour angle
  constants$RAD_PER_DAY <- 0.017214    # radians of Earth orbit per julian day
  constants$RAD_PER_DEG <- 0.01745329  # radians per degree
  constants$MIN_DECL <- -0.4092797     # minimum declination (radians)
  constants$DAYS_OFF <- 11.25          # julian day offset of winter solstice
  # Note:  Make sure that 3600 % SW_RAD_DT == 0
  constants$SW_RAD_DT <- 30.0          # timestep for radiation routine (seconds)
  
  constants$MA <- 28.9644e-3         # (kg mol-1) molecular weight of air
  constants$R <- 8.3143              # (m3 Pa mol-1 K-1) gas law constant
  constants$R_DRY <- 287             # (J / degC * kg) Gas constant of dry air
  constants$G_STD <- 9.80665         # (m s-2) standard gravitational accel.
  constants$P_STD <- 101325.0        # (Pa) standard pressure at 0. m elevation
  constants$T_STD <- 288.15          # (K) standard temp at 0. m elevation
  constants$CP <- 1010.0             # (J kg-1 K-1) specific heat of air
  
  constants$KELVIN <- 273.15
  constants$EPS <- 0.62196351
  constants$DAYS_PER_YEAR <- 365.25
  constants$SEC_PER_DAY <- 86400
  constants$SEC_PER_MIN <- 60
  constants$MIN_PER_HOUR <- 60
  constants$HOURS_PER_DAY <- 24
  constants$MIN_PER_DAY <- 1440
  constants$SEC_PER_HOUR <- 3600
  constants$STEFAN_B <- 5.669e-8      # (W m^-22 K^-4) Stefan Boltzmann constant
  
  constants$M_PER_KM <- 1000.0
  constants$MBAR_PER_BAR <- 1000.0
  constants$MM_PER_CM <- 10.0
  constants$MAX_PERCENT <- 100.0
  
  # parameters for the radiation algorithm
  # (dim) stands for dimensionless values
  constants$C <- 1.5             # (dim) radiation parameter
  constants$B0 <- 0.031          # (dim) radiation parameter
  constants$B1 <- 0.201          # (dim) radiation parameter
  constants$B2 <- 0.185          # (dim) radiation parameter
  constants$TBASE <- 0.870       # (dim) max inst. trans. 0m nadir dry atm
  constants$ABASE <- -6.1e-5     # (1/Pa) vapor pressure effect on transmittance
  constants$SC_INT <- 1.32       # (MJ/m2/day) snow correction intercept
  constants$SC_SLOPE <- 0.096    # (MJ/m2/day/cm) snow correction slope
  
  constants$EPS <- 0.62196351 # Ratio of molecular weights: M_water_vapor/M_dry_air
    
  # Allocate the radiation arrays
  tiny_step_per_day <- constants$SEC_PER_DAY / constants$SW_RAD_DT
  
  return(constants)
}
