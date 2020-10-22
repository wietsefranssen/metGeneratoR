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
  double **relhum_c = (double**)malloc(nx * sizeof(double));
  for (ix = 0; ix < nx; ix++) {
    relhum_c[ix] = (double*)malloc(ny * sizeof(double));
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

// [[Rcpp::export]]
NumericVector set_max_min_lonlat_cr(NumericVector tmin_map, NumericVector tmax_map, int yday, int nrec, NumericVector xybox) {
  float reslon = 0.5;
  float lon;
  int nx, ix;
  float slon = xybox[0];
  float elon = xybox[1];
  float slat = xybox[2];
  float elat = xybox[3];
  float reslat = 0.5;
  float lat;
  int ny, iy;
  // int nt;
  int nTinyStepsPerDay;
  
  
  int irec;
  
  ny = ((elat - slat) / reslat) + 1;
  nx = ((elon - slon) / reslon) + 1;
  
  NumericVector tair_map_r(nrec * ny * nx);
  IntegerVector dims(3);
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nrec;
  tair_map_r.attr("dim") = dims;
  
  // double tmin_hour_ix;
  // double tmax_hour_ix;
  // nt = nx;
  nTinyStepsPerDay = 360 / reslon; // 360 degress in lon direction
  
  // double *radfrac_c = (double*)malloc(nx * sizeof(double));
  double tmin_hour = -999;
  double tmax_hour = -999;;
  
  // Define and allocate
  double **rad_fract_map_org = (double**)malloc(ny * sizeof(double));
  for (iy = 0; iy < ny; iy++) {
    rad_fract_map_org[iy] = (double*)malloc(nTinyStepsPerDay * sizeof(double));
  }
  
  // Define and allocate
  double ***tair_map = (double***)malloc(nx * sizeof(double));
  for (ix = 0; ix < nx; ix++) {
    tair_map[ix] = (double**)malloc(ny * sizeof(double));
    for (iy = 0; iy < ny; iy++) {
      tair_map[ix][iy] = (double*)malloc(nrec * sizeof(double));
    }
  }
  
  double *Tair = (double*)malloc(nrec * sizeof(double));
  
  // run the function
  rad_fract_lats_c(rad_fract_map_org, nTinyStepsPerDay, yday, slat, elat);
  
  int ix_offset;
  float tmin_hour_new;
  float tmax_hour_new;
  for (iy = 0; iy < ny; iy++) {
    set_min_max_hour_c(rad_fract_map_org[iy], &tmin_hour, &tmax_hour, nTinyStepsPerDay);
    lat = slat + (iy*reslat);
    if (lat < 0 && tmin_hour == 99) {
      tmin_hour = 0;
      tmax_hour = 16;
    } 
    if (lat >= 0 && tmin_hour == 99) {
      tmin_hour = 11 - (24/(360/reslat));
      tmax_hour = 12 + (24/(360/reslat));
    }
    float iXoffset = (slon - -179.75) / reslon;
    
    // printf("iy: %d, tmin: %f, tmax: %f, radfrac: %f\n", iy, tmin_hour, tmax_hour, sum);
    for (ix = 0; ix < nx; ix++) {
      lon = slon + (ix*reslon);
      ix_offset = (lon / 0.5) - 0.5;
      // float iXoffset = (slon - -179.75) / reslon;
      tmin_hour_new = tmin_hour - (24.0/(float)nTinyStepsPerDay * (float)ix_offset);
      // printf("ix_offset: %d, iXoffset: %f\n", ix_offset,iXoffset);
      // if (tmin_hour_new<0) tmin_hour_new +=24;
      tmax_hour_new = tmax_hour - (24.0/(float)nTinyStepsPerDay * (float)ix_offset);
      // if (tmax_hour_new<0) tmax_hour_new +=24;
      // if (iy == 200)
      // printf("blaat: lon: %5.2f, iy: %d, ix: %d, ix_offset: %d, tmin: %f, tmax: %f, tminhour: %f\n", lon, iy, ix, ix_offset,tmin_hour_new, tmax_hour_new, tmin_hour);
      
      // double Tmin = 15;
      // double tmax = 30;
      // HourlyT_c(nrec, tmin_hour_new, Tmin, tmax_hour_new, tmax, Tair);
      HourlyT_c(nrec, tmin_hour_new, tmin_map[iy*nx+ix], tmax_hour_new, tmax_map[iy*nx+ix], Tair);
      
      // Copy back to R array
      for (irec = 0; irec < nrec; irec++) {
        tair_map[ix][iy][irec] = Tair[irec];
        // tair_map[ix][iy][irec] = tmin_map[iy*nx+ix];
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
  
  for (iy = 0; iy < ny; iy++) {
    free(rad_fract_map_org[iy]);
  }
  free(Tair);
  free(rad_fract_map_org);
  
  return tair_map_r;
}

//' @export
// [[Rcpp::export]]
NumericVector rad_map_final_cr(int nrec, int yday, double gmt_float, NumericVector xybox, NumericVector lats, NumericVector lons, float gmt_offset) {
  // Define and allocate
  int nx, ix;
  int ny, iy;
  int it;
  int irec;
  int idx;
  int iGmtOffset;
  int nTinyStepsPerDay;
  double dt = 30;
  double gmt_float_tmp;
  size_t count;
  double SEC_PER_DAY = 86400;
  ix = 0;
  float lon, lat;
  int irecit;
  int it_tmp;
  int tss;
  
  nTinyStepsPerDay = SEC_PER_DAY / dt;
  
  // dt = nTinyStepsPerDay / nrec;
  nx = (xybox[1] - xybox[0]) + 1;
  ny = (xybox[3] - xybox[2]) + 1;
  
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
  
  // Zero
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        rad_fract_map[ix][iy][irec] = 0;
      }
    }
  }
  
  // Define and allocate
  double *rad_fract = (double*)malloc(nTinyStepsPerDay * sizeof(double));
  for (it = 0; it < nTinyStepsPerDay; it++) {
    rad_fract[it] = 0;
  }
  
  
  // Calc gmt_offset
  if (gmt_offset < 0) gmt_offset += 24;
  gmt_offset += - (24/nrec) * 0.5;
  if (gmt_offset < 0) gmt_offset += 24;
  if (gmt_offset >= 24) gmt_offset -= 24;
  
  // run the function
  for (ix = 0; ix < nx; ix++) {
    iy = 0;
    lon = lons[nx*iy + ix];
    for (iy = 0; iy < ny; iy++) {
      lat = lats[nx*iy + ix];
      solar_geom_new_c(rad_fract, lat, yday, dt);
      it_tmp = floor(nTinyStepsPerDay * ( (lon + 180) / 360));
      if (it_tmp > nTinyStepsPerDay) it_tmp = it_tmp - nTinyStepsPerDay;
      it = floor(it_tmp + ( (nTinyStepsPerDay / 24) * gmt_offset));
      if (it > nTinyStepsPerDay) it = it - nTinyStepsPerDay;
      
      for (tss = 0; tss < nTinyStepsPerDay; tss++) {
        irec = floor(tss / (nTinyStepsPerDay / nrec));
        // irecit = floor(tss + it + (tss / nrec));
        irecit = floor(tss + it + ( (tss / nTinyStepsPerDay) * (nTinyStepsPerDay / nrec)));
        if (irecit >= nTinyStepsPerDay) {
          irecit = irecit - nTinyStepsPerDay;
        }
        rad_fract_map[ix][iy][irec] += rad_fract[irecit] * nrec;
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
  free(rad_fract);
  return rad_fract_map_r;
}

// [[Rcpp::export]]
NumericVector rad_map_lats_cr(int nt, int yday) {
  // Define and allocate
  float slat = -89.75;
  float elat = 89.75;
  float reslat = 0.5;
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
  rad_fract_lats_c(rad_fract_map, nt, yday, slat, elat);
  
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
