#include <Rcpp.h>
#include "header.h"

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

#include "Rcpp.h"        // R memory io
#include "Rdefines.h"        // R memory io
#include "Rmath.h"    // R math functions

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector main2(int yday) {
  int ntimesteps_per_day = 24;
  int nx = 720;
  int ny = 360;
  int nOutStepDay = 4;
  
  NumericVector final_result(nx * ny * nOutStepDay);
  
  float lat, lon;
  int ix, iy, it, ity;
  int irec;
  
  int hour_offset_int;
  int idx;
  size_t i;
  
  double *map_rad_tmp = (double*)malloc(ntimesteps_per_day * ny * sizeof(double));
  double *result = (double*)malloc(ntimesteps_per_day * sizeof(double));
  double *radfract_latlon_out = (double*)malloc(nx * ny * nOutStepDay * sizeof(double));
  
  for (i = 0; i < (nx * ny * nOutStepDay); i++) radfract_latlon_out[i] = 0;
  
  // map_rad_tmp
  for (iy = 0; iy < ny; iy++) {
    lat = -89.75 + ((iy)*0.5);
    solar_geom_c(result, lat, yday, ntimesteps_per_day);
    for (it = 0; it < ntimesteps_per_day; it++) {
      ity = (iy * ntimesteps_per_day) + it;
      map_rad_tmp[ity] = result[it];
      // printf("hoi: %d, %d, %d, %f, %f\n", ity, iy, it, map_rad_tmp[ity], result[it]);
    }
  }
  
  size_t id_maprad,iyrec,id_out,itt,ntt;
  // int nrOffsetSteps = 24;
  // 
  // ntt = ntimesteps_per_day / nOutStepDay;
  // 
  // for (ix = 0; ix < nx; ix++) {
  //   lon = (ix * 0.5) + -179.75;
  //   hour_offset_int =  floor((float)ix * ( (float)nrOffsetSteps / (float)nx) );
  //   // printf("hour_offset_int: %f, %d,  %d\n", lon, ix, hour_offset_int);
  //   for (iy = 0; iy < ny; iy++) {
  //     for (irec = 0; irec < nOutStepDay; irec++) {
  //       for (itt = 0; itt < ntt; itt++) {
  //         idx = (irec * ntt) + itt + hour_offset_int;
  //         if (idx >= ntimesteps_per_day) idx -= ntimesteps_per_day;
  //         int  itimesteps = (irec * ntt) + itt;
  //         // id_maprad = (iy * ntimesteps_per_day) + it;
  //         // id_maprad = (iy * ntimesteps_per_day) + it;
  //         id_maprad = (iy*ntimesteps_per_day) + itimesteps;
  //         
  //         iyrec = (iy * nOutStepDay) + irec;
  //         id_out = (ix * ny * nOutStepDay) + iyrec;
  //         
  //         radfract_latlon_out[id_out] += map_rad_tmp[id_maprad];// * nOutStepDay;
  //         printf("id_out: %zu, id_maprad: %zu, ix: %d, iy: %d, irec: %d, itt: %zu, itimesteps: %d, offset: %d, %f, %f\n", id_out, id_maprad, ix, iy, irec, itt, itimesteps, hour_offset_int, map_rad_tmp[id_maprad], radfract_latlon_out[id_out]);
  //       }
  //     }
  //   }
  // }
  
  
  // map_rad_tmp
  for (iy = 0; iy < ny; iy++) {
    lat = -89.75 + ((iy)*0.5);
    solar_geom_c(result, lat, yday, ntimesteps_per_day);
    for (it = 0; it < ntimesteps_per_day; it++) {
      ity = (iy * ntimesteps_per_day) + it;
      map_rad_tmp[ity] = result[it];
      printf("ity: %d, iy: %d, it: %d, %f, %f\n", ity, iy, it, map_rad_tmp[ity], result[it]);
    }
  }
  
  // map_rad_tmp
  // for (iy = 0; iy < ny; iy++) {
  //   lat = -89.75 + ((iy)*0.5);
  //   for (ix = 0; ix < nx; ix++) {
  //     it = floor((float)ix / ( (float)nx / (float)ntimesteps_per_day) );
  //     id_out = (iy * nx) + ix;
  //     id_maprad = (iy * ntimesteps_per_day) + it;
  //     radfract_latlon_out[id_out] = map_rad_tmp[id_maprad];
  //     printf("id_out: %d, id_maprad: %d, iy: %d, ix: %d, it: %d, %f, %f\n", id_out, id_maprad, iy, ix, it, map_rad_tmp[id_maprad],result[it]);
  //     // printf("ix: %d, nx,: %d, ntimesteps_per_day: %d, it: %d\n", ix, nx, ntimesteps_per_day,it);
  //   }
  // }
  
  
  // radfract_latlon(radfract_latlon_out, map_rad_tmp, nx, ny, ntimesteps_per_day, nOutStepDay);
  for (i = 0; i < (nx * ny * nOutStepDay); i++) final_result[i] = 0;
  for (i = 0; i < (nx * ny * nOutStepDay); i++) final_result[i] = radfract_latlon_out[i];
  // for (i = 0; i < (ntimesteps_per_day * ny); i++) final_result[i] = map_rad_tmp[i];
  
  free(radfract_latlon_out);
  free(map_rad_tmp);
  free(result);
  
  return final_result;
}

// [[Rcpp::export]]
NumericVector rad_map_final_cr(int nrec, int yday) {
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
  
  nt = 360;
  nt = 24;
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
  int hour_offset_int;
  int nrOffsetSteps = nt;
  int dt = nt/nrec;
  int idx = 0;
  size_t count = 0;
  for (ix = 0; ix < nx; ix++) {
    lon = (ix * 0.5) + -179.75;
    // hour_offset_int = ceil((ix + 1) * (nrOffsetSteps / nx)); // hour_offset<-0
    hour_offset_int =  floor((float)ix * ( (float)nrOffsetSteps / (float)nx) );
    // printf("offset: %d\n", hour_offset_int);
    for (iy = 0; iy < ny; iy++) {
        idx = hour_offset_int;
        for (irec = 0; irec < nrec; irec++) {
          for (int id = 0; id < dt; id++) {
            if (idx >= 24) idx = idx - nt;
            // rad_fract_map_r[count] = count;
            size_t qq = iy*nrec + irec;
            size_t qqq = ix*ny + qq;
            // rad_fract_map_r[qqq] += rad_fract_map_org[iy][idx] / nrec;
            rad_fract_map_r[count] += rad_fract_map_org[iy][idx] / nrec;
            idx++;  
            
          }
          count++;
          
        //     size_t qq = iy*nt + it; 
        // size_t qqq = ix*ny + qq; 
        // printf("nr: %zu, %zu\n",qq,qqq);
        // arr<-1:nrOffsetSteps
          
        // rad_fract_map_r[iy*nt + it] = rad_fract_map[iy][it];
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

int rad_map_final_c(double ***rad_fract_map_final, double **rad_fract_map, int nt, int yday) {
  
  float slat = -89.75;
  float elat = 89.75;
  float reslat = 0.5;
  float lat;
  int ny, iy;
  int it;
  
  ny = ((elat - slat) / reslat) + 1;
  
  // Define and allocate
  double *rad_fract = (double*)malloc(nt * sizeof(double));
  
  for (iy = 0; iy < ny; iy++) {
    lat = slat + ((iy)*0.5);
    
    // run the function
    // solar_geom_c(rad_fract, lat, yday, nt);
    for (it = 0; it < nt; it++) {
      // rad_fract_map[iy][it] = rad_fract[it];
      // *(rad_fract_map + iy*nt + it) = rad_fract[it];
    }
    
  } 
  // Free
  free(rad_fract);
  
  return 0;
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

/*** R
bb<-rad_map_lats_cr(24, 1)
*/

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

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */