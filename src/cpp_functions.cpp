#include <iostream>
#include <math.h>
#include <stdio.h>
#include "header.h"
#include <cmath>

float deg2rad(float deg) {
  return (deg * M_PI) / (180);
}

float potential_radiation(int hour, int minute, int yday, float lon, float lat, float timezone) {
  float terrain_slope=0;
  float terrain_slope_azimuth=0;
  float cloud_fraction=0;

  float solar_constant = 1367.0;
  float days_per_year = 365.25;
  float tropic_of_cancer = deg2rad(23.43697);
  // float tropic_of_cancer = 0.4090523;
  float solstice = 173.0;
  float pi = M_PI;
  
  int dates_hour = hour;
  int dates_minute = minute;
  int day_of_year = yday;
  
  // compute solar decline in rad
  float solar_decline = tropic_of_cancer * cos(2.0 * pi * (day_of_year - solstice) / days_per_year);
  // compute the sun hour angle in rad
  float standard_meridian = timezone * 15.;
  float delta_lat_time = (lon - standard_meridian) * 24. / 360.;
  float hour_angle = pi * (((dates_hour + dates_minute / 60. + delta_lat_time) / 12.) - 1.);
  
  // get solar zenith angle
  float cos_solar_zenith = (sin(solar_decline) * sin(deg2rad(lat))
                              + cos(solar_decline) * cos(deg2rad(lat)) * cos(hour_angle));
  if (cos_solar_zenith < 0) { cos_solar_zenith = 0; }
  float solar_zenith_angle = acos(cos_solar_zenith);
  
  // compute transmissivities for direct and diffus radiation using cloud fraction
  float transmissivity_direct = (0.6 + 0.2 * cos_solar_zenith) * (1.0 - cloud_fraction);
  float transmissivity_diffuse = (0.3 + 0.1 * cos_solar_zenith) * cloud_fraction;
  
  // modify solar constant for eccentricity
  float beta = 2. * pi * (day_of_year / days_per_year);
  float radius_ratio = (1.00011 + 0.034221 * cos(beta) + 0.00128 * sin(beta)
                           + 0.000719 * cos(2. * beta) + 0.000077 * sin(2 * beta));
  float solar_constant_times_radius_ratio = solar_constant * radius_ratio;
  
  float mu_tmp = cos(solar_decline) * sin(hour_angle) / sin(solar_zenith_angle);
  if (mu_tmp < -1.0) mu_tmp = -1.0;
  if (mu_tmp > 1.0) mu_tmp = 1.0;
  float mu = asin(mu_tmp);
  // float mu = asin(cos(solar_decline) * sin(hour_angle) / sin(solar_zenith_angle));
  float cosi = (cos(terrain_slope) * cos_solar_zenith
                   + sin(terrain_slope) * sin(solar_zenith_angle) * cos(mu - terrain_slope_azimuth));
  
  // get total shortwave radiation
  float direct_radiation = solar_constant_times_radius_ratio * transmissivity_direct * cosi;
  float diffuse_radiation = solar_constant_times_radius_ratio * transmissivity_diffuse * cos_solar_zenith;
  if (direct_radiation < 0) { direct_radiation = 0; }
  // printf("hoi%f\n",cos(solar_decline) );
  // printf("%f\n",cos(solar_decline) );
  // printf("%f\n",sin(hour_angle));
  // printf("%f\n",sin(solar_zenith_angle));
  // printf("%f\n",cos(solar_decline) * sin(hour_angle) / sin(solar_zenith_angle)  );
  // printf("mu: %f\n",asin(cos(solar_decline) * sin(hour_angle) / sin(solar_zenith_angle)));
  // 
  // printf("%f\n",mu_tmp);
  // printf("mu2: %f\n",(float)mu);
  // printf("%f\n",asin(-1.00000));
  // 
  return direct_radiation;
}

int set_min_max_hour_c(double *radfrac, double *tmin_hour, double *tmax_hour, int nx) {
  int ix;
  // int iy;
  int ix_prev;
  int tmin_ix = -999;
  int tmax_ix = -999;
  
  for (ix = 0; ix < (nx/2); ix++) {
    ix_prev = ix - 1;
    if (ix_prev < 0) { ix_prev = nx - 1;}
    if (radfrac[ix] > 0 && radfrac[ix_prev] <= 0)
    {
      tmin_ix = ix_prev;
      break;
    }
  }
  for (ix = (nx/2); ix < nx; ix++) {
    ix_prev = ix - 1;
    if (ix_prev < 0) { ix_prev = nx - 1;}
    if (radfrac[ix] <= 0 && radfrac[ix_prev] > 0)
    {
      tmax_ix = ix + 1;
      break;
    }
  }
  
  // convert ix to hour
  double partsPerHour = nx/24;
  *tmin_hour = -999;
  *tmax_hour = -999;
  
  if (tmin_ix >= 1 && tmax_ix >= 1) {
    *tmin_hour = tmin_ix / partsPerHour;
    *tmax_hour = tmax_ix / partsPerHour;
  }
  if (*tmin_hour >= 0 && *tmax_hour >= 0) {
    *tmax_hour = 0.67 * (*tmax_hour - *tmin_hour) + *tmin_hour;
  } else {
    /* arbitrarily set the min and max times to 2am and 2pm */
    // *tmin_hour = 2;
    // *tmax_hour = 14;
    *tmin_hour = 99;
    *tmax_hour = 99;
  }
  
  // printf("tmin_hour: %f tmax_hour: %f\n", *tmin_hour, *tmax_hour);
  return 0;
}

int solar_geom_c(double *rad_fract_per_timestep, float lat, int yday, int timesteps_per_day) {
  
  double dt = 30;
  double PII = 3.141593;
  double MIN_DECL = -0.4092797;
  double RAD_PER_DEG = 0.01745329;
  double SEC_PER_RAD = 13750.99;
  double SEC_PER_DAY = 86400;
  double DAYS_OFF = 11.25;
  double RAD_PER_DAY = 0.017214;
  double cosegeom, sinegeom, coshss, hss, daylength;
  double dir_beam_topa;
  
  int tinystepsperday, i, j, tss;
  double coslat, sinlat, dh, decl, cosdecl, sindecl;
  // double sum_trans = 0;
  double sum_flat_potrad = 0;
  double cosh, h, am, cza, dir_flat_topa;
  int tinystep;
  int istep;
  int ami;
  
  tinystepsperday = SEC_PER_DAY / dt;
  
  /* Allocate the radiation arrays */
  // double tiny_rad_fract_1day[tinystepsperday];
  double *tiny_rad_fract_1day = (double*)malloc(tinystepsperday * sizeof(double));
  
  for (i = 0; i < timesteps_per_day; i++) rad_fract_per_timestep[i] = 0;
  for (i = 0; i < tinystepsperday; i++) tiny_rad_fract_1day[i] = 0;
  
  /* optical airmass by degrees */
  double optam[21] = {2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07, 4.37, 4.72, 5.12, 5.60,
                      6.18, 6.88, 7.77, 8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00};
  
  /* precalculate the transcendentals */
  /* check for (+/-) 90 degrees latitude, throws off daylength calc */
  lat *= RAD_PER_DEG;
  if (lat > 1.570796)
    lat = 1.570796;
  if (lat < -1.570796)
    lat = -1.570796;
  coslat = cos(lat);
  sinlat = sin(lat);
  
  /*  Sub-daily time step and angular step */
  dh = dt / SEC_PER_RAD;
  
  /* Declination and quantities of interest */
  decl = MIN_DECL * cos(((yday - 1) + DAYS_OFF) * RAD_PER_DAY);
  cosdecl = cos(decl);
  sindecl = sin(decl);
  
  /* calculate daylength as a function of lat and decl */
  cosegeom = coslat * cosdecl;
  sinegeom = sinlat * sindecl;
  coshss = -(sinegeom) / cosegeom;
  if (coshss < -1.0)
    coshss = -1.0; /* 24-hr daylight */
  if (coshss > 1.0)
    coshss = 1.0; /* 0-hr daylight */
  hss = acos(coshss); /* hour angle at sunset (radians) */
  /* daylength (seconds) */
  daylength = 2.0 * hss * SEC_PER_RAD;
  if (daylength > 86400) daylength = 86400;
  
  /* extraterrestrial radiation perpendicular to beam, total over
   the timestep (J) */
  dir_beam_topa = (1368.0 + 45.5 * sin((2.0 * PII * (double) (yday - 1) / 365.25) + 1.7)) * dt;
  
  /* Set up angular calculations */
  for (h = -hss; h < hss; h += dh) {
    /* Cosine of the hour angle and solar zenith angle */
    cosh = cos(h);
    /* calculate cosine of solar zenith angle */
    cza = cosegeom * cosh + sinegeom;
    
    if (cza > 0.0) {
      /* When sun is above flat horizon do flat-surface
       calculations to determine daily total transmittance
       and save potential radiation for calculation of
       diffuse portion */
      
      /* potential radiation for this time period, flat surface,
       top of atmosphere */
      dir_flat_topa = dir_beam_topa * cza;
      /* determine optical air mass */
      am = 1.0 / (cza + 0.0000001);
      if (am > 2.9) {
        ami = (int) (acos(cza) / RAD_PER_DEG) - 69;
        if (ami < 0)
          ami = 0;
        if (ami > 20)
          ami = 20;
        am = optam[ami];
      }
      
      /* keep track of total potential radiation on a flat
       surface for ideal horizons */
      sum_flat_potrad += dir_flat_topa;
      
    }/* end if sun above ideal horizon */
      else dir_flat_topa = -1;
      
      /* start vic_change */
      tinystep = (12L * 3600L + h * SEC_PER_RAD) / dt;
      if (tinystep < 0)
        tinystep = 0;
      if (tinystep > tinystepsperday - 1)
        tinystep = tinystepsperday - 1;
      
      if (dir_flat_topa > 0)
        tiny_rad_fract_1day[tinystep] = dir_flat_topa;
      else
        tiny_rad_fract_1day[tinystep] = 0;
  }
  
  if (daylength && sum_flat_potrad > 0) {
    for (j = 0; j < tinystepsperday; j++) {
      tiny_rad_fract_1day[j] /= sum_flat_potrad;
    }
  }
  
  for (tss = 0; tss < tinystepsperday; tss++) {
    istep = floor(tss / (tinystepsperday / timesteps_per_day));
    rad_fract_per_timestep[istep] += tiny_rad_fract_1day[tss];
  }
  
  free(tiny_rad_fract_1day);
  
  return 0;
}

int rad_fract_lats_c(double **rad_fract_map, int nt, int yday, float slat, float elat) {
  
  // float slat = -89.75;
  // float elat = 89.75;
  float reslat = 0.25;
  float lat;
  int ny, iy;
  int it;
  
  ny = ((elat - slat) / reslat) + 1;
  
  // Define and allocate
  double *rad_fract = (double*)malloc(nt * sizeof(double));
  
  for (iy = 0; iy < ny; iy++) {
    lat = slat + ((iy)*0.25);
    
    // run the function
    solar_geom_c(rad_fract, lat, yday, nt);
    for (it = 0; it < nt; it++) {
      rad_fract_map[iy][it] = rad_fract[it];
      // *(rad_fract_map + iy*nt + it) = rad_fract[it];
    }
    
  } 
  // Free
  free(rad_fract);
  
  return 0;
}
