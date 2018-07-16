#include <math.h>
#include "Rcpp.h"        // R memory io
#include "Rdefines.h"        // R memory io

/*** SubRoutine Prototypes ***/
// double * solar_geom_c(float, int , int);
int solar_geom_c(double*, float, int , int);
int radfract_latlon(double *map_rad_tmp, double *radfrac, int nx, int ny, int nt, int nOutStepDay);
  
// // [[Rcpp::export]]
// int ccc();

// // [[Rcpp::export]]
// NumericVector main2();

// // [[Rcpp::export]]
// NumericVector timesTwo(NumericVector);