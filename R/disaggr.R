set_min_max_hour <- function(tiny_rad_fract, yday, n_days,
                             ts, params) {

  # def set_min_max_hour(tiny_rad_fract: np.array, yday: np.array, n_days: int,
  #                      ts: int, params: dict) -> Tuple[np.array]:
  #     """
  #     Determine the time at which min and max temp
  #     is reached for each day.
  #
  #     Parameters
  #     ----------
  #     tiny_rad_fract:
  #         Array of fraction of shortwave radiation received
  #         at a shortened timestep. This should be calculated
  #         by `metsim.physics.solar_geom`.
  #     yday:
  #         Array of day of year for each simulated day.
  #     n_days:
  #         Number of days in run.
  #     ts:
  #         Timestep of run.
  #     params:
  #         Dictionary of parameters to use.  Must contain
  #         `utc_offset` and `tmax_daylength_fraction`.
  #
  #     Returns
  #     -------
  #     (t_t_min, t_t_max):
  #         A tuple containing 2 timeseries, corresponding
  #         to time of min and max temp, respectively
  #     """
  load(file = "./hpc/outRadFractions_2880.Rdata")
  ## Select a day
  day <- 1
  tiny_rad_fract <- outRadFractions[200,day,]
  cnst <- metGen$cnst
  # calculate minute of sunrise and sunset for each day of the year
  rad_mask = 1 * (tiny_rad_fract > 0)
  # mask = np.diff(rad_mask)
  # rise_times = np.where(mask > 0)[1] * (cnst.SW_RAD_DT / cnst.SEC_PER_MIN)
  rise_times = np.where(mask > 0)[1] * (cnst$SW_RAD_DT / cnst$SEC_PER_MIN)
  # set_times = np.where(mask < 0)[1] * (cnst.SW_RAD_DT / cnst.SEC_PER_MIN)
  set_times = np.where(mask < 0)[1] * (cnst$SW_RAD_DT / cnst$SEC_PER_MIN)
  #
  # if params['utc_offset']:
  #     # not used elsewhere:
  #     # rad_fract_per_day = int(cnst.SEC_PER_DAY/cnst.SW_RAD_DT)
  #     utc_offset = int(((params.get("lon", 0) - params.get("theta_s", 0)) /
  #                       cnst.DEG_PER_REV) * cnst.MIN_PER_DAY)
  #     rise_times -= utc_offset
  #     set_times -= utc_offset
  #
  # # map the daily sunrise and sunset to a monotonic timseries (in minutes)
  # offset = np.arange(0, n_days * cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY,
  #                    cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY)
  # rise_times = rise_times[yday] + offset
  # set_times = set_times[yday] + offset
  #
  # # time of maximum and minimum temperature calculated thusly
  # t_t_max = (params['tmax_daylength_fraction'] * (set_times - rise_times) +
  #            rise_times) + ts
  # t_t_min = rise_times
  # return t_t_min, t_t_max
}