#include "header.h"
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector set_min_max_hour_cr(NumericVector radfrac, int nx) {
  NumericVector result(2);
  int ix;

  double *radfrac_c = (double*)malloc(nx * sizeof(double));
  double tmin_hour = -999;
  double tmax_hour = -999;;
  
  for (ix = 0; ix < nx; ix++) {
    radfrac_c[ix] = radfrac[ix];
  }
  
  set_min_max_hour_c(radfrac_c, &tmin_hour, &tmax_hour, nx);
    
  result[0] = tmin_hour;
  result[1] = tmax_hour;
  return result;
}
// qq<-set_max_min_hour_cr(solar_geom_cr(-89.75,1,720), 720)
// plot(solar_geom_cr(-89.75,1,720))
// plot(HourlyT(TmaxHour = qq[2],Tmax = 30,TminHour = 0,Tmin = 15)) 


// [[Rcpp::export]]
NumericVector set_max_min_lonlat_cr(int nx, int ny, int yday, int nrec) {
  NumericVector result(2);
  NumericVector tair_map_r(nrec * ny * nx);
  IntegerVector dims(3);
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nrec;
  tair_map_r.attr("dim") = dims;
  
  int ix, iy, nt;
  
  nt = nx;
  double *radfrac_c = (double*)malloc(nx * sizeof(double));
  double tmin_hour = -999;
  double tmax_hour = -999;;
  
  // Define and allocate
  double **rad_fract_map_org = (double**)malloc(ny * sizeof(double));
  for (iy = 0; iy < ny; iy++) {
    rad_fract_map_org[iy] = (double*)malloc(nt * sizeof(double));
  }
  
  // Define and allocate
  double ***tair_map = (double***)malloc(nx * sizeof(double));
  for (ix = 0; ix < nx; ix++) {
    tair_map[ix] = (double**)malloc(ny * sizeof(double));
    for (iy = 0; iy < ny; iy++) {
      tair_map[ix][iy] = (double*)malloc(nrec * sizeof(double));
    }
  }
  
  // run the function
  rad_fract_lats_c(rad_fract_map_org, nt, yday);
  
  iy=0;
  set_min_max_hour_c(rad_fract_map_org[iy], &tmin_hour, &tmax_hour, nx);

  double Tmin = 15;
  double Tmax = 30;
  // double TminHour = 8.1;
  // double TmaxHour = 14.1;
  double TminHour = 3;
  double TmaxHour = 12;
  
  double *Tair = (double*)malloc(nrec * sizeof(double));
  
  HourlyT_c(nrec, TminHour, Tmin, TmaxHour, Tmax, Tair);

  int irec;
  ix=0;
  iy=0;
  // for (ix = 0; ix < nx; ix++) {
    // for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        tair_map[ix][iy][irec] = Tair[irec];
        // printf("ttair: %f\n", Tair[irec]);
      }
    // }
  // }
  
  size_t count = 0;
  for (irec = 0; irec < nrec; irec++) {
    for (iy = 0; iy < ny; iy++) {
      for (ix = 0; ix < nx; ix++) {
        tair_map_r[count] = tair_map[ix][iy][irec];
        count++;
      }
    }
  }
  
  
  // Free
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      free(tair_map[ix][iy]);
    }
    free(tair_map[ix]);
  }
  free(tair_map);
  
  for (iy = 0; iy < ny; iy++) {
    free(rad_fract_map_org[iy]);
  }
  free(Tair);
  
  free(rad_fract_map_org);
  result[0] = tmin_hour;
  result[1] = tmax_hour;
  return tair_map_r;
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
