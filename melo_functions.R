# def potential_radiation(dates, lon, lat, timezone, terrain_slope=0, terrain_slope_azimuth=0,
#                         cloud_fraction=0, split=False):
#   """
# Calculate potential shortwave radiation for a specific location and time.
# 
# This routine calculates global radiation as described in:
# Liston, G. E. and Elder, K. (2006): A Meteorological Distribution System for
# High-Resolution Terrestrial Modeling (MicroMet), J. Hydrometeorol., 7, 217–234.
# 
# Corrections for eccentricity are carried out following:
# Paltridge, G.W., Platt, C.M.R., 1976. Radiative processes in Meteorology and Climatology.
# Elsevier Scientific Publishing Company, Amsterdam, Oxford, New York.
# 
# Parameters
# ----------
# dates : DatetimeIndex or array-like
# The dates for which potential radiation shall be calculated
# lon : float
# Longitude (degrees)
# lat : float
# Latitude (degrees)
# timezone : float
# Time zone
# terrain_slope : float, default 0
# Terrain slope as defined in Liston & Elder (2006) (eq. 12)
# terrain_slope_azimuth : float, default 0
# Terrain slope azimuth as defined in Liston & Elder (2006) (eq. 13)
# cloud_fraction : float, default 0
# Cloud fraction between 0 and 1
# split : boolean, default False
# If True, return a DataFrame containing direct and diffuse radiation,
# otherwise return a Series containing total radiation
# """

# rm(list = ls())
## Based on : https://journals.ametsoc.org/doi/pdf/10.1175/JHM486.1
## https://www.geosci-model-dev.net/9/2315/2016/gmd-2016-51.pdf
library(units)
library(lubridate)

deg2rad <- function(deg) {(deg * pi) / (180)}

potential_radiation <- function(hour,minute,yday,lon=8.86,lat= 51.00,timezone=1) {
  # potential_radiation <- function(dates='2014-01-01 12:00:00',lon=8.86,lat= 51.00,timezone=1) {
    # dates='2014-01-01 12:00:00'
  # lon=8.86
  # lat= 51.00
  # timezone=1
  terrain_slope=0
  terrain_slope_azimuth=0
  cloud_fraction=0 
  split=F
  
  solar_constant = 1367.
  days_per_year = 365.25
  tropic_of_cancer = deg2rad(23.43697)
  solstice = 173.0
  
  # dates = pd.DatetimeIndex(dates)
  # dates_hour = hour(dates)
  # dates_minute = minute(dates)
  # day_of_year = yday(dates)
  
  dates_hour = hour
  dates_minute = minute
  day_of_year = yday
  
  # compute solar decline in rad
  solar_decline = tropic_of_cancer * cos(2.0 * pi * (day_of_year - solstice) / days_per_year)
  
  # compute the sun hour angle in rad
  standard_meridian = timezone * 15.
  delta_lat_time = (lon - standard_meridian) * 24. / 360.
  hour_angle = pi * (((dates_hour + dates_minute / 60. + delta_lat_time) / 12.) - 1.)
  
  # get solar zenith angle
  cos_solar_zenith = (sin(solar_decline) * sin(deg2rad(lat))
                      + cos(solar_decline) * cos(deg2rad(lat)) * cos(hour_angle))
  cos_solar_zenith[cos_solar_zenith<0]<-0
  solar_zenith_angle = acos(cos_solar_zenith)
  
  # compute transmissivities for direct and diffus radiation using cloud fraction
  transmissivity_direct = (0.6 + 0.2 * cos_solar_zenith) * (1.0 - cloud_fraction)
  transmissivity_diffuse = (0.3 + 0.1 * cos_solar_zenith) * cloud_fraction
  
  # modify solar constant for eccentricity
  beta = 2. * pi * (day_of_year / days_per_year)
  radius_ratio = (1.00011 + 0.034221 * cos(beta) + 0.00128 * sin(beta)
                  + 0.000719 * cos(2. * beta) + 0.000077 * sin(2 * beta))
  solar_constant_times_radius_ratio = solar_constant * radius_ratio
  
  mu = asin(cos(solar_decline) * sin(hour_angle) / sin(solar_zenith_angle))
  cosi = (cos(terrain_slope) * cos_solar_zenith
          + sin(terrain_slope) * sin(solar_zenith_angle) * cos(mu - terrain_slope_azimuth))
  
  # get total shortwave radiation
  direct_radiation = solar_constant_times_radius_ratio * transmissivity_direct * cosi
  diffuse_radiation = solar_constant_times_radius_ratio * transmissivity_diffuse * cos_solar_zenith
  direct_radiation[direct_radiation<0]<-0
  
  return(direct_radiation)
}

##########################

disaggregate_radiation <- function(radiation=10.6625,hour,minute, yday,lon=8.86,lat= 51.00,timezone=1) {
  # disaggregate_radiation <- function(radiation=10.6625,dates='2014-01-01 12:00:00',lon=8.86,lat= 51.00,timezone=1) {
    # pot_rad <- potential_radiation(dates=dates,lon=lon,lat= lat,timezone=timezone)
  # pot_rad <- potential_radiation(hour = hour,minute = minute, yday =yday,lon=lon,lat= lat,timezone=timezone)
  pot_rad<-array(NA, dim = c(length(hour)))
  for (i in 1:length(hour)) {
    # hour = hour(dates[i])
    # minute = minute(dates[i])
    # yday =yday(dates[i])
    h <- hour[i]
    m <- minute[i]
    y <- yday[i]
    radtmp <- potential_radiation_cr(h, m, y, lon, lat, timezone)
    # radtmp <- potential_radiation_cr(hour = hour[i],minute = minute[i], yday =yday[i], lon=lon,lat= lat,timezone=timezone)
    # radtmp <- potential_radiation(hour = hour[i],minute = minute[i], yday =yday[i],lon=lon,lat= lat,timezone=timezone)
    pot_rad[i] <-radtmp
  }
  globalrad = radiation
  
  pot_rad_daily = mean(pot_rad)
  glob_disagg = pot_rad / pot_rad_daily * globalrad
  glob_disagg[glob_disagg < 1e-2] = 0.
  glob_disagg[is.na(glob_disagg)] = 0.
  
  # print(glob_disagg)
  return(glob_disagg)
}
