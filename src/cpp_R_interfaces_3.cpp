#include "header.h"
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rad_map_final3_cr(int nrec, int yday) {
  // Define and allocate
  float slon = -179.75;
  float elon = 179.75;
  float reslon = 0.5;
  float lon;
  int nx, ix;
  float slat = -89.75;
  float elat = 89.75;
  float reslat = 0.5;
  float res = reslat;
  float lat;
  int ny, iy;
  int it;
  size_t count = 0;
  
  int irec;
  int nTinyStepsPerDay;
  
  ny = ((elat - slat) / reslat) + 1;
  nx = ((elon - slon) / reslon) + 1;
  nTinyStepsPerDay = 360 / res;
  
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
    rad_fract_map_org[iy] = (double*)malloc(nTinyStepsPerDay * sizeof(double));
  }
  
  for (iy = 0; iy < ny; iy++) {
    for (it = 0; it < nTinyStepsPerDay; it++) {
      rad_fract_map_org[iy][it] = 0;
    }
  }

  // run the function
  rad_fract_lats_c(rad_fract_map_org, nTinyStepsPerDay, yday);
  
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        rad_fract_map[ix][iy][irec] = 0;
      }
    }
  }
  
  double gmt_float = 11;
  double gmt_float_tmp;
  int iGmtOffset, iLonOffset;
  
  // ## Check and correct gmt_float 
  // if(gmt_float < -12 || gmt_float > 12) stop("cannot be lower than -12 and higher than 12")
  if (gmt_float < 0) gmt_float = gmt_float + nTinyStepsPerDay;
  
  // ## Define the index offset based on the gmt offset
  gmt_float_tmp = ( gmt_float + (24/nrec)/2 ); // + ((nx/2)/nrec);
  iGmtOffset = gmt_float_tmp * (nTinyStepsPerDay/24);
  iGmtOffset = 270;
  // printf("gmt: %d\n", iGmtOffset);
  
  
  // // OLD *******************************************************************
  int offset;
  int dt = nTinyStepsPerDay/nrec;
  int idx = 0;
  int idx_tmp = 0;
  int gmt = 0;
  count = 0;
  for (ix = 0; ix < nx; ix++) {
    idx_tmp = ix;
    if (idx_tmp >= nx) idx_tmp = idx_tmp - nx;
    idx_tmp = idx_tmp + iGmtOffset ;
    if (idx_tmp >= nx) idx_tmp = idx_tmp - nx;
    idx = ix;
    // printf("ix: %d, idx: %d\n", ix, idx_tmp);
    
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        idx = (irec * (nx/nrec)) + idx_tmp;
        // idx = (irec * (nx/nrec)) + ix;
        if (idx >= nx) idx = idx - nx;
        rad_fract_map[ix][iy][irec] = rad_fract_map_org[iy][idx] * nx;
      }
    }
  }
  // // OLD *******************************************************************
  
  

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
