# # data_daily =
# method='sine_min_max'
# min_max_time='fix'
# mod_nighttime=FALSE
# # max_delta=None
# # mean_course=None
# # sun_times=None
# 
# hours_per_day = 24
# method = 'sine_min_max'
# min_max_time = 'fix'
# 
# data_daily.temp <- 273.
# data_daily.tmin <- 271.05
# data_daily.tmax <- 276.45


temp_disagg = numeric(hours_per_day)

if (method == 'sine_min_max') {
  # for this option assume time of minimum and maximum and fit cosine function through minimum and maximum temperatures
  default_shift_hours = 2
  
  daylength_thres = 3
  # min / max hour during polar night assumption
  min_loc_polar = 6
  max_loc_polar = 18
  
  locdf <- NULL
  locdf$index  <- NULL ## DATE
  locdf$min_loc  <- NULL
  locdf$max_loc  <- NULL
  locdf$min_val_before  <- NULL
  locdf$min_val_cur  <- NULL
  locdf$min_val_next  <- NULL
  locdf$max_val_before  <- NULL
  locdf$max_val_cur  <- NULL
  locdf$max_val_next  <- NULL
  locdf$mean_val_cur  <- NULL
  
  if (min_max_time == 'fix') {
    # take fixed location for minimum and maximum
    locdf$min_loc = 7
    locdf$max_loc = 14
  } else if (min_max_time == 'sun_loc') {
    #   # take location for minimum and maximum by sunrise / sunnoon + 2h
    #   locdf.min_loc = sun_times.sunrise.round()  # sun rise round to full hour
    #   locdf.max_loc = sun_times.sunnoon.round() + default_shift_hours  # sun noon round to full hour + fix 2h
  } else if (min_max_time == 'sun_loc_shift') {
    #   # take location for minimum and maximum by sunrise / sunnoon + monthly delta
    #   locdf.min_loc = sun_times.sunrise.round()  # sun rise round to full hour
    #   locdf.max_loc = (sun_times.sunnoon + max_delta[locdf.index.month].values).round()  # sun noon + shift derived from observed hourly data, round to full hour
    # 
    #   pos = locdf.min_loc > locdf.max_loc
    #   locdf.loc[pos, 'max_loc'] = sun_times.sunnoon[pos].round() + default_shift_hours  # standard shift in this case
  }
  
  locdf$min_val_cur = data_daily.tmin
  locdf$max_val_cur = data_daily.tmax
  locdf$mean_val_cur = data_daily.temp
  locdf$min_val_next = data_daily.tmin
  locdf$max_val_next = data_daily.tmax
  # locdf.loc[locdf.index[-1], 'min_val_next'] = locdf.min_val_cur.iloc[-1]
  # locdf.loc[locdf.index[-1], 'max_val_next'] = locdf.max_val_cur.iloc[-1]
  locdf$min_val_before = data_daily.tmin
  locdf$max_val_before = data_daily.tmax
  # locdf.loc[locdf.index[0], 'min_val_before'] = locdf.min_val_cur.iloc[0]
  # locdf.loc[locdf.index[0], 'max_val_before'] = locdf.max_val_cur.iloc[0]
  
  locdf_day = locdf
  # locdf = locdf.reindex(temp_disagg.index, method='ffill')
  
  print(locdf)
  # # whenever we are before the maximum for the current day, use minimum value of current day for cosine function fitting
  # # once we have passed the maximum value use the minimum for next day to ensure smooth transitions
  min_val = locdf.min_val_next
  # min_val[min_val.index.hour < locdf.max_loc] = locdf.min_val_cur
  # 
  # # whenever we are before the minimum for the current day, use maximum value of day before for cosine function fitting
  # # once we have passed the minimum value use the maximum for the current day to ensure smooth transitions
  max_val = locdf.max_val_cur
  # max_val[max_val.index.hour < locdf.min_loc] = locdf.max_val_before
  # 
  # temp_disagg = pd.Series(index=min_val.index)
  # 
  if (method == 'sine_min_max') {
    delta_val = max_val - min_val
    v_trans = min_val + delta_val / 2.
    
    if(mod_nighttime) {
      #   before_min = locdf.index.hour <= locdf.min_loc
      # between_min_max = (locdf.index.hour > locdf.min_loc) & (locdf.index.hour < locdf.max_loc)
      # after_max = locdf.index.hour >= locdf.max_loc
      # temp_disagg[before_min]      = v_trans + delta_val / 2. * np.cos(np.pi / (hours_per_day - (locdf.max_loc - locdf.min_loc)) * (hours_per_day - locdf.max_loc + locdf.index.hour))
      # temp_disagg[between_min_max] = v_trans + delta_val / 2. * np.cos(1.25 * np.pi + 0.75 * np.pi / (locdf.max_loc - locdf.min_loc) * (locdf.index.hour - locdf.min_loc))
      # temp_disagg[after_max]       = v_trans + delta_val / 2. * np.cos(np.pi / (hours_per_day - (locdf.max_loc - locdf.min_loc)) * (locdf.index.hour - locdf.max_loc))
    } else {
      # temp_disagg = v_trans + (delta_val / 2.) * cos(2 * pi / hours_per_day * (locdf$index.hour - locdf$max_loc))
      for (i in 1:hours_per_day) 
        temp_disagg[i] = v_trans + (delta_val / 2.) * cos(2 * pi / hours_per_day * ((i-1) - locdf$max_loc))
    }
  } else if (method == 'sine_mean') {
     # dtr = locdf$max_val_cur - locdf$min_val_cur
    # temp_disagg[:] = locdf.mean_val_cur + dtr / 2. * np.cos(2 * np.pi / hours_per_day * (locdf.index.hour - locdf.max_loc))
  }
  print(temp_disagg)
  plot(c(0:23),temp_disagg)
  # polars = sun_times.daylength < daylength_thres
  # if polars.sum() > 0:
  #   # during polar night, no diurnal variation of temperature is applied
  #   # instead the daily average calculated using tmin and tmax is applied
  #   polars_index_hourly = melodist.util.hourly_index(polars[polars].index)
  # temp_disagg.loc[polars_index_hourly] = np.nan
  # 
  # avg_before = (locdf_day.min_val_before + locdf_day.max_val_before) / 2.
  # avg_cur = (locdf_day.min_val_cur + locdf_day.max_val_cur) / 2.
  # getting_warmers = polars &  (avg_before <= avg_cur)
  # getting_colders = polars & ~(avg_before <= avg_cur)
  # 
  # getting_warmers_min_loc = pd.DatetimeIndex([ts.replace(hour=min_loc_polar) for ts in getting_warmers[getting_warmers].index])
  # getting_warmers_max_loc = pd.DatetimeIndex([ts.replace(hour=max_loc_polar) for ts in getting_warmers[getting_warmers].index])
  # temp_disagg[getting_warmers_min_loc] = locdf_day.min_val_cur[getting_warmers].values
  # temp_disagg[getting_warmers_max_loc] = locdf_day.max_val_cur[getting_warmers].values
  # 
  # getting_colders_min_loc = pd.DatetimeIndex([ts.replace(hour=min_loc_polar) for ts in getting_colders[getting_colders].index])
  # getting_colders_max_loc = pd.DatetimeIndex([ts.replace(hour=max_loc_polar) for ts in getting_colders[getting_colders].index])
  # temp_disagg[getting_colders_min_loc] = locdf_day.max_val_cur[getting_colders].values
  # temp_disagg[getting_colders_max_loc] = locdf_day.min_val_cur[getting_colders].values
  # 
  # temp_polars = temp_disagg.loc[polars_index_hourly].copy()
  # transition_days = polars[polars.diff() == True].astype(int) # 0 where transition from polar to "normal" mode, 1 where transition from normal to polar
  # 
  # if len(transition_days) > 0:
  #   polar_to_normal_days = transition_days.index[transition_days == 0]
  # normal_to_polar_days = transition_days.index[transition_days == 1] - pd.Timedelta(days=1)
  # add_days = polar_to_normal_days.union(normal_to_polar_days)
  # 
  # temp_polars = temp_polars.append(temp_disagg[melodist.util.hourly_index(add_days)]).sort_index()
  # 
  # for day in polar_to_normal_days:
  #   min_loc = int(locdf.loc[day].min_loc)
  # temp_polars[day.replace(hour=0):day.replace(hour=min_loc) - pd.Timedelta(hours=1)] = np.nan
  # temp_polars[day.replace(hour=min_loc)] = locdf.min_val_cur[day]
  # 
  # for day in normal_to_polar_days:
  #   max_loc = int(locdf.loc[day].max_loc)
  # temp_polars[day.replace(hour=max_loc) + pd.Timedelta(hours=1):day.replace(hour=23)] = np.nan
  # 
  # temp_interp = temp_polars.interpolate(method='linear', limit=23)
  # temp_disagg[temp_interp.index] = temp_interp
}
# elif method == 'mean_course_min_max':
#   data_daily_as_hourly = data_daily.reindex(temp_disagg.index, method='ffill', limit=23)
# 
# df = pd.DataFrame(index=temp_disagg.index)
# df['normval'] = mean_course.unstack().loc[list(zip(df.index.month, df.index.hour))].values
# df['tmin'] = data_daily_as_hourly.tmin
# df['tmax'] = data_daily_as_hourly.tmax
# 
# temp_disagg[:] = df.normval * (df.tmax - df.tmin) + df.tmin
# elif method == 'mean_course_mean':
#   data_daily_as_hourly = data_daily.reindex(temp_disagg.index, method='ffill', limit=23)
# dtr = data_daily_as_hourly.tmax - data_daily_as_hourly.tmin
# mc = pd.Series(index=temp_disagg.index)
# mean_course_zeromean = mean_course - mean_course.mean()  # shift mean course so that the daily mean is 0
# mc[:] = mean_course_zeromean.unstack().loc[list(zip(temp_disagg.index.month, temp_disagg.index.hour))].values
# temp_disagg[:] = data_daily_as_hourly.temp + dtr * mc
# 
# return temp_disagg
# 
# 
# def get_shift_by_data(temp_hourly, lon, lat, time_zone):
#   '''function to get max temp shift (monthly) by hourly data
#     
#     Parameters
#     ----
#     hourly_data_obs : observed hourly data 
#     lat :             latitude in DezDeg
#     lon :             longitude in DezDeg
#     time_zone:        timezone
#     '''
# daily_index = temp_hourly.resample('D').mean().index
# sun_times = melodist.util.get_sun_times(daily_index, lon, lat, time_zone)
# 
# idxmax = temp_hourly.groupby(temp_hourly.index.date).idxmax()
# idxmax.index = pd.to_datetime(idxmax.index)
# max_temp_hour_obs = idxmax.dropna().apply(lambda d: d.hour)
# max_temp_hour_pot = sun_times.sunnoon
# max_delta = max_temp_hour_obs - max_temp_hour_pot
# mean_monthly_delta = max_delta.groupby(max_delta.index.month).mean()
# 
# return mean_monthly_delta
