#include <math.h>
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

/*** SubRoutine Prototypes ***/
int solar_geom_c(double*, float, int, double);
int rad_map_final_c(double ***rad_fract_map_final, double **rad_fract_map, int nt, int yday);
int set_min_max_hour_c(double *radfrac, double *tmin_hour, double *tmax_hour, int nx);

int set_sunrise_sunset_hour_c(double *radfract, double *sunrise, double *noon, double *sunset, int nx);
int set_tmin_tmax_hour_c(double sunrise, double noon, double sunset, double *tmin_hour, double *tmax_hour, int nx);

void set_max_min_hour(double *hourlyrad,  int tmaxhour, int tminhour); // weg?
void HourlyT_c(int nrec,
             double TmaxHour, 
             double Tmax, 
             double TminHour,
             double Tmin, 
             double *Tair);

double svp(double temp);
double sh2vp_c(double q, double t, double p);