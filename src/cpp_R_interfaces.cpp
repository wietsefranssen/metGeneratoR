#include "header.h"
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector set_vp_cr(NumericVector tair_r, NumericVector relhum_r, int nx, int ny, int nrec) {
  int irec,ix,iy;
  NumericVector vp_r(nx*ny*nrec);
  IntegerVector dims(3);
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nrec;
  vp_r.attr("dim") = dims;
  
  // Define and allocate
  float **relhum_c = (float**)malloc(nx * sizeof(float*));
  for (ix = 0; ix < nx; ix++) {
    relhum_c[ix] = (float*)malloc(ny * sizeof(float));
  }
  
  // Copy from R array
  size_t count = 0;
  for (iy = 0; iy < ny; iy++) {
    for (ix = 0; ix < nx; ix++) {
      relhum_c[ix][iy] = relhum_r[count];
      count++;
    }
  }
  
  // Do the thing
  count = 0;
  for (irec = 0; irec < nrec; irec++) {
    for (iy = 0; iy < ny; iy++) {
      for (ix = 0; ix < nx; ix++) {
        vp_r[count] = svp(tair_r[count]) * relhum_c[ix][iy] / 1000 ;
        count++;
      }
    }
  }
  
  // Free
  for (ix = 0; ix < nx; ix++) {
    free(relhum_c[ix]);
  }
  free(relhum_c);
  
  
  
  return vp_r; 
}  

// [[Rcpp::export]]
NumericVector sh2vp(NumericVector q, NumericVector p) {
  size_t i;
  size_t n = q.length();
  
  NumericVector t = p;
  
  if(q.length() != t.length() || q.length() != p.length()) {
    // if(q.length() != t.length() || q.length() != p.length()) {
    printf("q and p should have the same length!\n");
    // printf("q, t and p should all have the same length!\n");
    return(1);
  }
  
  NumericVector vp(n);
  vp.attr("dim") = q.attr("dim");
  
  for (i = 0; i < n; i++) {
    vp[i] = sh2vp_c(q[i], t[i], p[i]);
  }
  
  return vp;
}  

//' @export
// [[Rcpp::export]]
NumericVector calc_tas_cr(NumericVector tmin_hour_map, NumericVector tmax_hour_map, NumericVector tmin_map, NumericVector tmax_map, int yday, int nrec, NumericVector xybox) {
  float lon;
  int nx, ix;
  int ny, iy;
  float lat;
  int irec;
  
  nx = (xybox[1] - xybox[0]) + 1;
  ny = (xybox[3] - xybox[2]) + 1;
  
  NumericVector tair_map_r(nrec * ny * nx);
  IntegerVector dims(3);
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nrec;
  tair_map_r.attr("dim") = dims;
  
  // Define and allocate
  float ***tair_map = (float***)malloc(nx * sizeof(float**));
  for (ix = 0; ix < nx; ix++) {
    tair_map[ix] = (float**)malloc(ny * sizeof(float*));
    for (iy = 0; iy < ny; iy++) {
      tair_map[ix][iy] = (float*)malloc(nrec * sizeof(float));
    }
  }
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        tair_map[ix][iy][irec] = 0;
      }
    }
  }
  
  float *Tair = (float*)malloc(nrec * sizeof(float));
  
  
  for (iy = 0; iy < ny; iy++) {
    for (ix = 0; ix < nx; ix++) {
      if (!isnan(tmin_hour_map[iy*nx+ix])) {
        HourlyT_c(nrec, tmin_hour_map[iy*nx+ix], tmin_map[iy*nx+ix], tmax_hour_map[iy*nx+ix], tmax_map[iy*nx+ix], Tair);
      }
      for (irec = 0; irec < nrec; irec++) {
        tair_map[ix][iy][irec] = Tair[irec];
      }
    }
  }
  
  // Copy back to R array
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
  free(Tair);
  
  return tair_map_r;
}

//' @export
// [[Rcpp::export]]
NumericVector rad_map_final_cr(int nrec, int yday, NumericVector xybox, NumericVector lats, NumericVector lons, float gmt_offset) {
  // Define and allocate
  int nx, ix;
  int ny, iy;
  int it;
  int irec;
  int idx;
  int iGmtOffset;
  int nTinyStepsPerDay;
  float dt = 30;
  size_t count;
  int SEC_PER_DAY = 86400;
  ix = 0;
  float lon, lat;
  int irecit;
  int it_tmp;
  int tss;
  
  nTinyStepsPerDay = SEC_PER_DAY / dt;
  
  // dt = nTinyStepsPerDay / nrec;
  nx = (xybox[1] - xybox[0]) + 1;
  ny = (xybox[3] - xybox[2]) + 1;
  
  NumericVector rad_fract_map_r((nrec+2) * ny * nx);
  IntegerVector dims(3);
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nrec+2;
  rad_fract_map_r.attr("dim") = dims;
  
  // Define and allocate
  float ***rad_fract_map = (float***)malloc(nx * sizeof(float**));
  for (ix = 0; ix < nx; ix++) {
    rad_fract_map[ix] = (float**)malloc(ny * sizeof(float*));
    for (iy = 0; iy < ny; iy++) {
      rad_fract_map[ix][iy] = (float*)malloc(nrec * sizeof(float));
    }
  }
  float **t_min_hour_map = (float**)malloc(nx * sizeof(float*));
  for (ix = 0; ix < nx; ix++) {
    t_min_hour_map[ix] = (float*)malloc(ny * sizeof(float));
  }
  float **t_max_hour_map = (float**)malloc(nx * sizeof(float*));
  for (ix = 0; ix < nx; ix++) {
    t_max_hour_map[ix] = (float*)malloc(ny * sizeof(float));
  }
  
  // Zero
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      t_min_hour_map[ix][iy] = 0;
      t_max_hour_map[ix][iy] = 0;
      for (irec = 0; irec < nrec; irec++) {
        rad_fract_map[ix][iy][irec] = 0;
      }
    }
  }
  
  // Define and allocate
  float *rad_fract = (float*)malloc(nTinyStepsPerDay * sizeof(float));
  for (it = 0; it < nTinyStepsPerDay; it++) {
    rad_fract[it] = 0;
  }
  
  // Calc gmt_offset
  gmt_offset += - (24/nrec) * 0.5;
  gmt_offset += 12;
  gmt_offset -= floor(gmt_offset / 24) * 24;
  
  // run the function
  for (ix = 0; ix < nx; ix++) {
    iy = 0;
    for (iy = 0; iy < ny; iy++) {
      lon = lons[nx*iy + ix];
      lat = lats[nx*iy + ix];
      solar_geom_c(rad_fract, lat, yday, dt);
      it_tmp = floor(nTinyStepsPerDay * ( (lon + 180) / 360));
      it = floor(it_tmp + ( (nTinyStepsPerDay / 24) * (gmt_offset)));
      it -= floor(it / nTinyStepsPerDay) * nTinyStepsPerDay;
      
      for (tss = 0; tss < nTinyStepsPerDay; tss++) {
        irec = floor(tss / (nTinyStepsPerDay / nrec));
        irecit = floor(tss + it + ( (tss / nTinyStepsPerDay) * (nTinyStepsPerDay / nrec)));
        irecit -= floor(irecit / nTinyStepsPerDay) * nTinyStepsPerDay;
        
        rad_fract_map[ix][iy][irec] += rad_fract[irecit] * nrec;
      }
      
      float sunrise, noon, sunset;
      float tmin_hour = -999;
      float tmax_hour = -999;
      
      set_sunrise_sunset_hour_c(rad_fract_map[ix][iy], &sunrise, &noon, &sunset, nrec);
      if (ix == 1 && iy == 1) printf("iy: %i, ix: %i, sunrise: %f, noon: %f, sunset: %f\n",iy, ix,sunrise, noon,sunset);
      
      int allZero = 1;
      for (int i = 0; i < 24; i++) {
        if (rad_fract_map[ix][iy][i] > 0) allZero = 0;
      }
      
      if (allZero) {
        sunrise = (it / (nTinyStepsPerDay/24)) + 6;
        noon =  (it / (nTinyStepsPerDay/24)) + 12;
        sunset =  (it / (nTinyStepsPerDay/24)) + 18;
        sunrise -= floor(sunrise / 24) * 24;
        noon -= floor(noon / 24) * 24;
        sunset -= floor(sunset / 24) * 24;
      }
      
      set_tmin_tmax_hour_c(sunrise, noon, sunset, &t_min_hour_map[ix][iy], &t_max_hour_map[ix][iy], nrec);
      // set_tmin_tmax_hour_c(sunrise, noon, sunset, &tmin_hour, &tmax_hour, nrec);
      // if (ix == 1 && iy == 1) printf("iy: %i, ix: %i, tmin_hour: %f, tmax_hour: %f\n",iy, ix,tmin_hour,tmax_hour);
      if (ix == 1 && iy == 1) printf("iy: %i, ix: %i, tmin_hour: %f, tmax_hour: %f\n",iy, ix,t_min_hour_map[ix][iy],t_max_hour_map[ix][iy]);
      if (allZero) {
        t_min_hour_map[ix][iy] = -(it / (nTinyStepsPerDay/24)) + (12-1);
        t_min_hour_map[ix][iy] -= floor(t_min_hour_map[ix][iy] / 24) * 24;
        t_max_hour_map[ix][iy] =  -(it / (nTinyStepsPerDay/24)) + (12);
        t_max_hour_map[ix][iy] -= floor(t_max_hour_map[ix][iy] / 24) * 24;
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
  // TminHour
  for (iy = 0; iy < ny; iy++) {
    for (ix = 0; ix < nx; ix++) {
      rad_fract_map_r[count] = t_min_hour_map[ix][iy];
      count++;
    }
  }
  // TmaxHour
  for (iy = 0; iy < ny; iy++) {
    for (ix = 0; ix < nx; ix++) {
      rad_fract_map_r[count] = t_max_hour_map[ix][iy];
      count++;
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
  
  for (ix = 0; ix < nx; ix++) {
    free(t_min_hour_map[ix]);
    free(t_max_hour_map[ix]);
  }
  free(t_min_hour_map);
  free(t_max_hour_map);
  
  free(rad_fract);
  return rad_fract_map_r;
}
