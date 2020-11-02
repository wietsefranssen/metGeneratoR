#include <math.h>
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

/*** SubRoutine Prototypes ***/
int solar_geom_c(float*, float, int, float);
int set_sunrise_sunset_hour_c(float *radfract, float *sunrise, float *noon, float *sunset, int nx);
int set_tmin_tmax_hour_c(float sunrise, float noon, float sunset, float *tmin_hour, float *tmax_hour, int nx);
void HourlyT_c(int nrec, float TmaxHour, float Tmax, float TminHour, float Tmin, float *Tair);
float svp(float temp);
float sh2vp_c(float q, float t, float p);
