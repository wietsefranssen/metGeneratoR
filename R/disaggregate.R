source("./R/constants.R")

shortwave <- function(sw_rad, daylength, day_of_year,
                      tiny_rad_fract, params) {


  #   """
  # Disaggregate shortwave radiation down to a subdaily timeseries.
  #
  # Parameters
  # ----------
  # sw_rad:
  # Daily incoming shortwave radiation
  # daylength:
  # List of daylength time for each day of year
  # day_of_year:
  # Timeseries of index of days since Jan-1
  # tiny_rad_fract:
  # Fraction of the daily potential radiation
  # during a radiation time step defined by SW_RAD_DT
  # params:
  # Dictionary of parameters from the MetSim object
  #
  # Returns
  # -------
  # disaggrad:
  # A sub-daily timeseries of shortwave radiation.
  # """
  ts <- params$time_step
  ts_hourly <- ts / cnst$MIN_PER_HOUR
  ts_per_day <- (cnst$HOURS_PER_DAY * cnst$MIN_PER_HOUR / ts)
  tmp_rad <- sw_rad * ts_per_day
  # tmp_rad <- (sw_rad * daylength) / (cnst$SEC_PER_HOUR * ts_hourly)

  n_days <- length(tmp_rad)
  disaggrad <- rep(0, n_days * ts_per_day)
  rad_fract_per_day <- (cnst$SEC_PER_DAY / cnst$SW_RAD_DT)
  # if params['utc_offset']:
  #         utc_offset = int(((params.get("lon", 0) - params.get("theta_s", 0)) /
  #                           cnst.DEG_PER_REV) * rad_fract_per_day)
  #         tiny_rad_fract = np.roll(tiny_rad_fract.flatten(), -utc_offset)
  #     else: {
  # tiny_rad_fract <- tiny_rad_fract.flatten()
  # }bincount
  chunk_size <- (ts * (cnst$SEC_PER_MIN / cnst$SW_RAD_DT))
  ts_id <- rep(c(1:ts_per_day),each=chunk_size)
  for (day in c(1:n_days)) {
    day<-2
    rad <- tiny_rad_fract[day, 1:rad_fract_per_day]
    st<-((day-1) * ts_per_day)+1
    end<-day * ts_per_day
    dslice <- c(st:end)
    rad_chunk <- rep(NA,ts_per_day)
    for (tss in 1:ts_per_day) {
      rad_chunk[tss] <- sum(rad[ts_id==tss])
    }
    disaggrad[dslice] <- rad_chunk * tmp_rad[day]
  }
  return (disaggrad)
}

shortwave_day <- function(sw_rad, daylength,
                          tiny_rad_fract, params) {


  #   """
  # Disaggregate shortwave radiation down to a subdaily timeseries.
  #
  # Parameters
  # ----------
  # sw_rad:
  # Daily incoming shortwave radiation
  # daylength:
  # List of daylength time for each day of year
  # tiny_rad_fract:
  # Fraction of the daily potential radiation
  # during a radiation time step defined by SW_RAD_DT
  # params:
  # Dictionary of parameters from the MetSim object
  #
  # Returns
  # -------
  # disaggrad:
  # A sub-daily timeseries of shortwave radiation.
  # """
  ts <- params$time_step
  ts_hourly <- ts / cnst$MIN_PER_HOUR
  ts_per_day <- (cnst$HOURS_PER_DAY * cnst$MIN_PER_HOUR / ts)
  tmp_rad <- sw_rad * ts_per_day
  # tmp_rad <- (sw_rad * daylength) / (cnst$SEC_PER_HOUR * ts_hourly)

  disaggrad <- rep(0, ts_per_day)
  rad_fract_per_day <- (cnst$SEC_PER_DAY / cnst$SW_RAD_DT)
  # if params['utc_offset']:
  #         utc_offset = int(((params.get("lon", 0) - params.get("theta_s", 0)) /
  #                           cnst.DEG_PER_REV) * rad_fract_per_day)
  #         tiny_rad_fract = np.roll(tiny_rad_fract.flatten(), -utc_offset)
  #     else: {
  # tiny_rad_fract <- tiny_rad_fract.flatten()
  # }bincount
  chunk_size <- (ts * (cnst$SEC_PER_MIN / cnst$SW_RAD_DT))
  ts_id <- rep(c(1:ts_per_day),each=chunk_size)

  rad <- tiny_rad_fract
  rad_chunk <- rep(NA,ts_per_day)
  for (tss in 1:ts_per_day) {
    rad_chunk[tss] <- sum(rad[ts_id==tss])
  }
  disaggrad <- rad_chunk * tmp_rad

  return (disaggrad)
}

tiny_rad_fract_chunk <- function(daylength, tiny_rad_fract, params) {
  ts <- params$time_step
  ts_hourly <- ts / cnst$MIN_PER_HOUR
  ts_per_day <- (cnst$HOURS_PER_DAY * cnst$MIN_PER_HOUR / ts)

  rad_fract_per_day <- (cnst$SEC_PER_DAY / cnst$SW_RAD_DT)
  chunk_size <- (ts * (cnst$SEC_PER_MIN / cnst$SW_RAD_DT))
  ts_id <- rep(c(1:ts_per_day),each=chunk_size)

  rad_chunk <- array(NA, dim=c(366,ts_per_day))
  for (day in c(1:366)) {
    rad <- tiny_rad_fract[day, 1:rad_fract_per_day]
    for (tss in 1:ts_per_day) {
      rad_chunk[day,tss] <- sum(rad[ts_id==tss])
    }
  }

  return (rad_chunk)
}

rad_fract_chunk_small <- function(dayLength, RadFractions, params) {
  
  ########## Calculate Rad_chunk
  ts <- params$time_step
  ts_per_day <- (cnst$HOURS_PER_DAY * cnst$MIN_PER_HOUR / ts)
  
  ny <- dim(outRadFractions)[1]
  
  rad_chunk <- array(0, dim = c(ny, 366, ts_per_day) )
  
  for (iy in 1:ny) {
    rad_chunk[iy,,] <- tiny_rad_fract_chunk(outDaylength[iy, ], outRadFractions[iy,,], params)
  }
  rad_chunk[is.na(rad_chunk)] <- 0
  # plot(rad_chunk[,1 ,3 ])
  return(rad_chunk)
}

shortwave_day_fast <- function(sw_rad, rad_chunk, daylength) {
  ts <- params$time_step
  ts_hourly <- ts / cnst$MIN_PER_HOUR
  ts_per_day <- (cnst$HOURS_PER_DAY * cnst$MIN_PER_HOUR / ts)

  tmp_rad <- sw_rad * ts_per_day
  # tmp_rad <- (sw_rad * daylength) / (cnst$SEC_PER_HOUR * ts_hourly)
  disaggrad <- rep(0, ts_per_day)
  disaggrad <- rad_chunk * tmp_rad

  return (disaggrad)
}
