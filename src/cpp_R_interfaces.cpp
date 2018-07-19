#include "header.h"
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector set_max_min_hour_cr(NumericVector hourlyrad_r, int nrec, int ix) {
  // /* optical airmass by degrees */
  // double optam[21] = {2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07, 4.37, 4.72, 5.12, 5.60,
  //                     6.18, 6.88, 7.77, 8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00};
  // 
  // int nrec = 24;
  int tmaxhour = -999; 
  int tminhour = -999;
  
  // Define and allocate
  NumericVector result(2);
  double *hourlyrad = (double*)malloc(24 * sizeof(double));
  
  int nt = nrec;
  int base_offset_gmt = 0;
  base_offset_gmt =  (nt / nrec) * (base_offset_gmt + 9);
  // int ix = 2;
  int nx = 720;
  int offset;
  offset = (floor((float)ix * ( (float)nt / (float)nx) )) + base_offset_gmt;
  if (offset >= nrec) offset = offset - nrec;
  printf("offset: %d\n", offset);
  
  
  
  
  
  // Pass array from R NumericVector 
  for (int i = 0; i < nrec; i++) hourlyrad[i] = hourlyrad_r[i];
  
  // run the function
  // set_max_min_hour(hourlyrad, tmaxhour, tminhour);
  int hour, prev_hour;
  for (hour = 0; hour < nrec; hour++) {
    prev_hour = hour - 1;
    if (prev_hour < 0) prev_hour = nrec - 1; 
    if (hourlyrad[hour] > 0 && hourlyrad[prev_hour] <= 0)
    {
      tminhour = prev_hour;
    }
    if (hourlyrad[hour] <= 0 && hourlyrad[prev_hour] > 0)
    {
      tmaxhour = hour + 1;
    }
  }
  if (tminhour == -999) {
    // Define lowest point
    double n = hourlyrad[0];
    for (hour = 0; hour < nrec; hour++) {
      if(hourlyrad[hour] < n) {
        n = hourlyrad[hour];
        tmaxhour = hour;
      }
    }
    tminhour = tmaxhour + 1;
    if (tminhour >= nrec) tminhour = 0;
  }
  // Moet dit?:
  if (tminhour >= 0 && tmaxhour >= 0) {
    tmaxhour = 0.67 * (tmaxhour - tminhour) + tminhour;
    
  } else {
    /* arbitrarily set the min and max times to 2am and 2pm */
    tminhour = 2;
    tmaxhour = 14;
  }
  
  // Pass array back to R NumericVector 
  result[0] = tminhour;
  result[1] = tmaxhour;
  
  // Free
  free(hourlyrad);
  return result;
}


// [[Rcpp::export]]
NumericVector rad_map_final_cr(int nrec, int yday, int nx_parts, double gmt_float) {
  // Define and allocate
  float slon = -179.75;
  float elon = 179.75;
  float reslon = 0.5;
  float lon;
  int nx, ix;
  float slat = -89.75;
  float elat = 89.75;
  float reslat = 0.5;
  float lat;
  int ny, iy;
  int nt, it;
  int ixx ;
  
  int irec;
  
  nt = nx_parts;
  ny = ((elat - slat) / reslat) + 1;
  nx = ((elon - slon) / reslon) + 1;
  
  NumericVector rad_fract_map_r(nrec * ny * nx);
  IntegerVector dims(3);
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nrec;
  rad_fract_map_r.attr("dim") = dims;
  
  // Define and allocate
  double ***rad_fract_map = (double***)malloc(nx * sizeof(double));
  for (ix = 0; ix < nx; ix++) {
    rad_fract_map[ix] = (double**)malloc(ny * sizeof(double));
    for (iy = 0; iy < ny; iy++) {
      rad_fract_map[ix][iy] = (double*)malloc(nrec * sizeof(double));
    }
  }
  
  // Define and allocate
  double **rad_fract_map_org = (double**)malloc(ny * sizeof(double));
  for (iy = 0; iy < ny; iy++) {
    rad_fract_map_org[iy] = (double*)malloc(nt * sizeof(double));
  }
  
  for (iy = 0; iy < ny; iy++) {
    for (it = 0; it < nt; it++) {
      rad_fract_map_org[iy][it] = 0;
    }
  }
  
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        rad_fract_map[ix][iy][irec] = 0;
      }
    }
  }
  
  // run the function
  rad_fract_lats_c(rad_fract_map_org, nt, yday);
  
  // Pass array back to R NumericVector
  // int base_offset = 1;
  int offset;
  int dt = nt/nrec;
  int idx = 0;
  size_t count = 0;
  int gmt = 0;
  int HoursPerDay = 24;
  
  int offset_index = (nt / HoursPerDay) * (gmt + 12);
  // int offset_index =  (nrec / HoursPerDay) * (gmt + 12);
  
  ///////////////
  // double gmt_float = 14;
  double gmt_float_tmp;
  int iGmtOffset, iLonOffset;
  int nTinyStepsPerDay =nt;
  // ## Check and correct gmt_float 
  // if(gmt_float < -12 || gmt_float > 12) stop("cannot be lower than -12 and higher than 12")
  if (gmt_float < 0) gmt_float = gmt_float + nTinyStepsPerDay;
  
  // ## Define the index offset based on the gmt offset
  // gmt_float_tmp = ( gmt_float + (24/nrec)/2 ); // + ((nx/2)/nrec);
  // iGmtOffset = gmt_float_tmp * (nTinyStepsPerDay/24);
  // iGmtOffset = 360;
  /////////////////
  gmt_float_tmp = gmt_float * (nTinyStepsPerDay/24);
  iGmtOffset = gmt_float_tmp + ((24/nrec) * -15 + 360);
  // printf("gmt: %d\n", iGmtOffset);
  if (iGmtOffset < 0) iGmtOffset = iGmtOffset + nTinyStepsPerDay;
  // printf("gmt: %d\n", iGmtOffset);
  
  
  // Do everything for rec 0
  for (ix = 0; ix < nx; ix++) {
    idx = (floor((float)ix * ( (float)nt / (float)nx) )) + iGmtOffset;
    for (int id = 0; id < dt; id++) {
      if (idx >= nt) idx -= nt;
      for (iy = 0; iy < ny; iy++) {
        rad_fract_map[ix][iy][0] += rad_fract_map_org[iy][idx] * nrec;
      }
      idx++;
    }
  }

  // Copy the map of rec 0 and move it for the other recs
  for (irec = 1; irec < nrec; irec++) {
    for (ix = 0; ix < nx; ix++) {
      ixx = ix - ( irec * (nx/nrec));
      if (ixx < 0) ixx += nx;
      for (iy = 0; iy < ny; iy++) {
        rad_fract_map[ixx][iy][irec] = rad_fract_map[ix][iy][0];
      }  
    }
  }
  
  count = 0;
  for (irec = 0; irec < nrec; irec++) {
    for (iy = 0; iy < ny; iy++) {
      for (ix = 0; ix < nx; ix++) {
        rad_fract_map_r[count] = rad_fract_map[ix][iy][irec];
        count++;
      }
    }
  }
  
  // Free
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      free(rad_fract_map[ix][iy]);
    }
    free(rad_fract_map[ix]);
  }
  free(rad_fract_map);
  
  for (iy = 0; iy < ny; iy++) {
    free(rad_fract_map_org[iy]);
  }
  free(rad_fract_map_org);
  return rad_fract_map_r;
}

// [[Rcpp::export]]
NumericVector rad_map_lats_cr(int nt, int yday) {
  // Define and allocate
  float slat = -89.75;
  float elat = 89.75;
  float reslat = 0.5;
  float lat;
  int ny, iy;
  int it;
  
  ny = ((elat - slat) / reslat) + 1;
  
  NumericVector rad_fract_map_r(nt * ny);
  IntegerVector dims(2);
  dims[0] = nt;
  dims[1] = ny;
  rad_fract_map_r.attr("dim") = dims;
  
  // Define and allocate
  double **rad_fract_map = (double**)malloc(ny * sizeof(double));
  for (iy = 0; iy < ny; iy++) {
    rad_fract_map[iy] = (double*)malloc(nt * sizeof(double));
  }
  
  // run the function
  rad_fract_lats_c(rad_fract_map, nt, yday);
  
  // Pass array back to R NumericVector
  for (iy = 0; iy < ny; iy++) {
    for (it = 0; it < nt; it++) {
      rad_fract_map_r[iy*nt + it] = rad_fract_map[iy][it];
    }
  }
  
  // Free
  for (iy = 0; iy < ny; iy++) {
    free(rad_fract_map[iy]);
  }
  free(rad_fract_map);
  
  return rad_fract_map_r;
}

// [[Rcpp::export]]
NumericVector solar_geom_cr(float lat, int yday, int timesteps_per_day) {
  
  // Define and allocate
  NumericVector result(timesteps_per_day);
  double *result_c = (double*)malloc(timesteps_per_day * sizeof(double));
  
  // run the function
  solar_geom_c(result_c, lat, yday, timesteps_per_day);
  
  // Pass array back to R NumericVector 
  for (int i = 0; i < timesteps_per_day; i++) result[i] = result_c[i];
  
  // Free
  free(result_c);
  return result;
}
