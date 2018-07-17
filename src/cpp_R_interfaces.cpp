#include "header.h"
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rad_map_final_cr(int nrec, int yday, int nx_parts) {
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
  
  int irec;
  
  nt = nx_parts;
  ny = ((elat - slat) / reslat) + 1;
  nx = ((elon - slon) / reslon) + 1;
  
  NumericVector rad_fract_map_r(nrec * ny * nx);
  IntegerVector dims(3);
  dims[0] = nrec;
  dims[1] = ny;
  dims[2] = nx;
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
  
  // run the function
  rad_fract_lats_c(rad_fract_map_org, nt, yday);
  
  // Pass array back to R NumericVector
  // int base_offset = 1;
  int offset;
  int nrOffsetSteps = nt;
  int dt = nt/nrec;
  int idx = 0;
  size_t count = 0;
  int base_offset_gmt = 0;
  
  base_offset_gmt =  (nt / 24) * (base_offset_gmt + 9);
  for (ix = 0; ix < nx; ix++) {
    offset =  (floor((float)ix * ( (float)nrOffsetSteps / (float)nx) )) + base_offset_gmt;
    idx = offset;
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        for (int id = 0; id < dt; id++) {
          if (idx >= nt) idx = idx - nt;
          rad_fract_map_r[count] += rad_fract_map_org[iy][idx] * nrec;
          idx++;  
        }
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
