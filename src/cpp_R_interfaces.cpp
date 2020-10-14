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
NumericVector rad_map_final_cr(int nrec, int yday, double gmt_float, NumericVector xybox) {
  // Define and allocate
  int nx, ix;
  int ny, iy;
  int it;
  // int ixx ;
  int irec;
  int idx;
  int iGmtOffset;
  int nTinyStepsPerDay;
  int dt;
  double gmt_float_tmp;
  size_t count;
  
  float reslon = 0.5;
  float reslat = 0.5;
  float slon = xybox[0];
  float elon = xybox[1];
  float slat = xybox[2];
  float elat = xybox[3];
  
  nTinyStepsPerDay = 360 / reslon; // 360 degress in lon direction
  dt = nTinyStepsPerDay / nrec;
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
    rad_fract_map_org[iy] = (double*)malloc(nTinyStepsPerDay * sizeof(double));
  }
  
  // Zero
  for (iy = 0; iy < ny; iy++) {
    for (it = 0; it < nTinyStepsPerDay; it++) {
      rad_fract_map_org[iy][it] = 0;
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
  
  // run the function
  rad_fract_lats_c(rad_fract_map_org, nTinyStepsPerDay, yday, slat, elat);
  
  // ## Check and correct gmt_float 
  // if(gmt_float < -12 || gmt_float > 12) stop("cannot be lower than -12 and higher than 12")
  if (gmt_float < 0) gmt_float = gmt_float + nTinyStepsPerDay;
  
  // ## Define the index offset based on the gmt offset
  gmt_float_tmp = gmt_float * (nTinyStepsPerDay/24);
  iGmtOffset = gmt_float_tmp + ((24/nrec) * -15 + 360);
  if (iGmtOffset < 0) iGmtOffset = iGmtOffset + nTinyStepsPerDay;
  // printf("gmt: %d\n", iGmtOffset);
  
  // iXoffset is for a smaller domain...
  float iXoffset = (slon - -179.75) / reslon;
  for (ix = 0; ix < nx; ix++) {
    // idx = (floor((float)ix * ( (float)nTinyStepsPerDay / (float)nx) )) + iGmtOffset;
    idx = floor((float)ix) + iGmtOffset + ((float)1 * (float)iXoffset);
    for (irec = 0; irec < nrec; irec++) {
      for (int id = 0; id < dt; id++) {
        if (idx >= nTinyStepsPerDay) idx -= nTinyStepsPerDay;
        for (iy = 0; iy < ny; iy++) {
          rad_fract_map[ix][iy][irec] += rad_fract_map_org[iy][idx] * nrec;
        }
        idx++;
      }
    }
  }
  
  // Do everything for rec 0
  // float iXoffset = (slon - -179.75) / reslon;
  // printf("iXoffset: %f\n", iXoffset);
  // // Do everything for rec 0
  // for (ix = 0; ix < nx; ix++) {
  //   // idx = (floor((float)ix * ( (float)nTinyStepsPerDay / (float)nx) )) + iGmtOffset;
  //   idx = floor((float)ix) + iGmtOffset + ((float)1 * (float)iXoffset);
  //   for (int id = 0; id < dt; id++) {
  //     if (idx >= nTinyStepsPerDay) idx -= nTinyStepsPerDay;
  //     for (iy = 0; iy < ny; iy++) {
  //       rad_fract_map[ix][iy][0] += rad_fract_map_org[iy][idx] * nrec;
  //     }
  //     idx++;
  //   }
  // }
  // 
  // // Copy the map of rec 0 and move it for the other recs
  // for (irec = 1; irec < nrec; irec++) {
  //   for (ix = 0; ix < nx; ix++) {
  //     ixx = ix - ( irec * (nx/nrec));
  //     if (ixx < 0) ixx += nx;
  //     for (iy = 0; iy < ny; iy++) {
  //       rad_fract_map[ixx][iy][irec] = rad_fract_map[ix][iy][0];
  //     }  
  //   }
  // }
  
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
//' @export
// [[Rcpp::export]]
NumericVector rad_map_final_2dll_cr(int nrec, int yday, double gmt_float, NumericVector xybox, NumericVector lats) {
  // Define and allocate
  int nx, ix;
  int ny, iy;
  int it;
  // int ixx ;
  int irec;
  int idx;
  int iGmtOffset;
  int nTinyStepsPerDay;
  int dt;
  double gmt_float_tmp;
  size_t count;
  // 
  // float reslon = 0.5;
  // float reslat = 0.5;
  // float slon = xybox[0];
  // float elon = xybox[1];
  // float slat = xybox[2];
  // float elat = xybox[3];
  // 
  // nTinyStepsPerDay = 360 / reslon; // 360 degress in lon direction
  // dt = nTinyStepsPerDay / nrec;
  // ny = ((elat - slat) / reslat) + 1;
  // nx = ((elon - slon) / reslon) + 1;

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
  
  // // Define and allocate
  // double **rad_fract_map_org = (double**)malloc(ny * sizeof(double));
  // for (iy = 0; iy < ny; iy++) {
  //   rad_fract_map_org[iy] = (double*)malloc(nTinyStepsPerDay * sizeof(double));
  // }
  // 
  // // Zero
  // for (iy = 0; iy < ny; iy++) {
  //   for (it = 0; it < nTinyStepsPerDay; it++) {
  //     rad_fract_map_org[iy][it] = 0;
  //   }
  // }
  
  // Zero
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      for (irec = 0; irec < nrec; irec++) {
        rad_fract_map[ix][iy][irec] = 0;
      }
    }
  }
  
  // // run the function new!
  // for (ix = 0; ix < nx; ix++) {
  //   for (iy = 0; iy < ny; iy++) {
  //     for (irec = 0; irec < nrec; irec++) {
  //       rad_fract_map[ix][iy][irec] = solar_geom_c(rad_fract, lat, yday, nt);
  //     }
  //   }
  // }
  // 
  // // // run the function
  // // rad_fract_lats_c(rad_fract_map_org, nTinyStepsPerDay, yday, slat, elat);
  // 
  // // ## Check and correct gmt_float 
  // // if(gmt_float < -12 || gmt_float > 12) stop("cannot be lower than -12 and higher than 12")
  // if (gmt_float < 0) gmt_float = gmt_float + nTinyStepsPerDay;
  // 
  // // ## Define the index offset based on the gmt offset
  // gmt_float_tmp = gmt_float * (nTinyStepsPerDay/24);
  // iGmtOffset = gmt_float_tmp + ((24/nrec) * -15 + 360);
  // if (iGmtOffset < 0) iGmtOffset = iGmtOffset + nTinyStepsPerDay;
  // 
  // // iXoffset is for a smaller domain...
  // float iXoffset = (slon - -179.75) / reslon;
  // for (ix = 0; ix < nx; ix++) {
  //   idx = floor((float)ix) + iGmtOffset + ((float)1 * (float)iXoffset);
  //   for (irec = 0; irec < nrec; irec++) {
  //     for (int id = 0; id < dt; id++) {
  //       if (idx >= nTinyStepsPerDay) idx -= nTinyStepsPerDay;
  //       for (iy = 0; iy < ny; iy++) {
  //         rad_fract_map[ix][iy][irec] += rad_fract_map_org[iy][idx] * nrec;
  //       }
  //       idx++;
  //     }
  //   }
  // }
  
  // Do everything for rec 0
  // float iXoffset = (slon - -179.75) / reslon;
  // printf("iXoffset: %f\n", iXoffset);
  // // Do everything for rec 0
  // for (ix = 0; ix < nx; ix++) {
  //   // idx = (floor((float)ix * ( (float)nTinyStepsPerDay / (float)nx) )) + iGmtOffset;
  //   idx = floor((float)ix) + iGmtOffset + ((float)1 * (float)iXoffset);
  //   for (int id = 0; id < dt; id++) {
  //     if (idx >= nTinyStepsPerDay) idx -= nTinyStepsPerDay;
  //     for (iy = 0; iy < ny; iy++) {
  //       rad_fract_map[ix][iy][0] += rad_fract_map_org[iy][idx] * nrec;
  //     }
  //     idx++;
  //   }
  // }
  // 
  // // Copy the map of rec 0 and move it for the other recs
  // for (irec = 1; irec < nrec; irec++) {
  //   for (ix = 0; ix < nx; ix++) {
  //     ixx = ix - ( irec * (nx/nrec));
  //     if (ixx < 0) ixx += nx;
  //     for (iy = 0; iy < ny; iy++) {
  //       rad_fract_map[ixx][iy][irec] = rad_fract_map[ix][iy][0];
  //     }  
  //   }
  // }
  
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
  
  // for (iy = 0; iy < ny; iy++) {
  //   free(rad_fract_map_org[iy]);
  // }
  // free(rad_fract_map_org);
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
