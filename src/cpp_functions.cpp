#include <iostream>
#include <math.h>
#include <stdio.h>
#include "header.h"

int set_sunrise_sunset_hour_c(float *radfract, float *sunrise, float *noon, float *sunset, int nx) {
  bool first;
  // sunrise
  first = 1;
  for (int i = 0; i < 24; i++) {
    if (i == 0) {
      if (radfract[i] > radfract[23]) {
        if (first) *sunrise = (float)i;
        first = 0;
      }
    } else {
      if (radfract[i] > radfract[i-1]) {
        if (first) {
          *sunrise = (float)i;
          break;
        }
      } else {
        if (first == 0) first = 1;
      }
    }
  }
  
  // Noon
  int imax = 0;
  for (int i = 0; i < 24; i++) {
    if (radfract[imax] < radfract[i]) {
      imax = i;
    }
  }
  *noon = (float)imax;
  
  // sunset
  for (int i = 23; i >= 0; i--) {
    if (i > 0) {
      if (radfract[i-1] > radfract[i]) {
        *sunset = (float)i;
        break;
      }
    } else {
      if (radfract[23] > radfract[i]) {
        *sunset = (float)i;
        break;
      }
    }
  }   

  return 0;
}
int set_tmin_tmax_hour_c(float sunrise, float noon, float sunset, float *tmin_hour, float *tmax_hour, int nx) {
  *tmin_hour = sunrise;
  *tmax_hour = (sunset - noon);
   if (*tmax_hour < 0) *tmax_hour += 24;
  
  *tmax_hour = noon + (*tmax_hour * 0.33);
  *tmax_hour -= floor(*tmax_hour / 24) * 24;
  if (*tmax_hour < 0) printf("tmax_hour: %f (%f, %f)\n", *tmax_hour, sunset, noon);
  // *tmax_hour = sunset;
  // } else {
  //   /* arbitrarily set the min and max times to 2am and 2pm */
  //   // */
  //   *tmin_hour = 2;
  //   *tmax_hour = 14;
  // }
  
  // printf("tmin_hour: %f tmax_hour: %f\n", *tmin_hour, *tmax_hour);
  return 0;
}

int solar_geom_c(float *tiny_rad_fract_1day, float lat, int yday, float dt = 30) {
  
  // float dt = 30;
  float PII = 3.141593;
  float MIN_DECL = -0.4092797;
  float RAD_PER_DEG = 0.01745329;
  float SEC_PER_RAD = 13750.99;
  float SEC_PER_DAY = 86400;
  float DAYS_OFF = 11.25;
  float RAD_PER_DAY = 0.017214;
  float cosegeom, sinegeom, coshss, hss, daylength;
  float dir_beam_topa;
  
  int tinystepsperday, i, j, tss;
  float coslat, sinlat, dh, decl, cosdecl, sindecl;
  // float sum_trans = 0;
  float sum_flat_potrad = 0;
  float cosh, h, am, cza, dir_flat_topa;
  int tinystep;
  int istep;
  int ami;
  
  tinystepsperday = SEC_PER_DAY / dt;
  
  for (i = 0; i < tinystepsperday; i++) tiny_rad_fract_1day[i] = 0;
  
  /* optical airmass by degrees */
  float optam[21] = {2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07, 4.37, 4.72, 5.12, 5.60,
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
  dir_beam_topa = (1368.0 + 45.5 * sin((2.0 * PII * (float) (yday - 1) / 365.25) + 1.7)) * dt;
  
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
  
  return 0;
}
