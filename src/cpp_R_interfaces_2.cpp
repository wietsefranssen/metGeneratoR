#include "header.h"
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rad_map_final2_cr(int nrec, int yday) {
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
  printf("gg %f %d %d \n", rad_fract_map_org[0][0], nTinyStepsPerDay, ny);
  
  // run the function
  rad_fract_lats_c(rad_fract_map_org, nTinyStepsPerDay, yday);

  // int base_offset = 1;
  // int offset;
  // int dt = nt/nrec;
  // int idx = 0;
  // size_t count = 0;
  // int gmt = 0;
  // int HoursPerDay = 24;
  
  double gmt_float = 0;
  double gmt_float_tmp;
  int iGmtOffset, iLonOffset;
  
  // ## Check and correct gmt_float 
  // if(gmt_float < -12 || gmt_float > 12) stop("cannot be lower than -12 and higher than 12")
  if (gmt_float < 0) gmt_float = gmt_float + nTinyStepsPerDay;
  
  // ## Define the index offset based on the gmt offset
  gmt_float_tmp = ( gmt_float + (24/nrec)/2 ) + ((nx/2)/nrec);
  iGmtOffset = gmt_float_tmp * (nTinyStepsPerDay/24);

  
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        rad_fract_map[ix][iy][irec] = 0;
      }
    }
  }

  // // OLD *******************************************************************
  int offset;
  int dt = nTinyStepsPerDay/nrec;
  int idx = 0;
  int gmt = 0;
  int HoursPerDay = 24;
  int offset_index = (nTinyStepsPerDay / HoursPerDay) * (gmt + 12);
  count = 0;
  for (ix = 0; ix < nx; ix++) {
    lon = (ix * res) + -179.75;
    iLonOffset = (( ((float)lon / (float)res) + (float)res ) )/ (720.0/(float)nTinyStepsPerDay);
    if (iLonOffset < 0) iLonOffset = iLonOffset + nTinyStepsPerDay;
    // idx = (floor((float)ix * ( (float)nTinyStepsPerDay / (float)nx) )) + offset_index + iLonOffset + iGmtOffset;
    // idx = (floor((float)ix * ( (float)nTinyStepsPerDay / (float)nx) )) + iLonOffset + iGmtOffset;
    idx = iLonOffset + iGmtOffset;
    
    // idx = (floor((float)ix * ( (float)nTinyStepsPerDay / (float)nx) ) + iGmtOffset);
    // idx = (floor((float)ix * ( (float)nt / (float)nx) )) + offset_index;
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        for (int id = 0; id < dt; id++) {
          if (idx >= nTinyStepsPerDay) idx = idx - nTinyStepsPerDay;
          rad_fract_map[ix][iy][irec] += rad_fract_map_org[iy][idx] * nrec;
          idx++;
        }
        count++;
      }
    }
  }
  // // OLD *******************************************************************
  
  
  // // // NEW *******************************************************************
  // double *radfrac_new = (double*)malloc(nTinyStepsPerDay * sizeof(double));
  // int ii = 0;
  // for (ix = 0; ix < nx; ix++) {
  //   lon = (ix * res) + -179.75;
  //   // ## Define the index offset longitude location
  //   iLonOffset = (( (lon / res) + res ) )/ (720/nTinyStepsPerDay);
  //   if (iLonOffset < 0) iLonOffset = iLonOffset + nTinyStepsPerDay;
  //   for (iy = 0; iy < ny; iy++) {
  //     for (int i = 0; i < nTinyStepsPerDay; i++) {
  //       int iNew = i + iGmtOffset + iLonOffset;
  //       if (iNew >= nTinyStepsPerDay) iNew = (iNew - nTinyStepsPerDay);
  //       radfrac_new[iNew] = rad_fract_map_org[iy][i];
  //       // printf("radfrac: %f\n",radfrac_new[iNew]);
  //     }
  // 
  //     for (int i = 0; i < (nTinyStepsPerDay); i++) {
  //       for (irec = 0; irec < nrec; irec++) {
  //         rad_fract_map[ix][iy][irec] = radfrac_new[i];
  //       }
  //     }
  //   }
  // }
  // // NEW *******************************************************************
  
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
