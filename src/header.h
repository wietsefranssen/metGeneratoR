#include <math.h>
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

/*** SubRoutine Prototypes ***/
int solar_geom_c(double*, float, int , int);
int solar_geom_new_c(double*, float, int, double);
int rad_fract_lats_c(double **rad_fract_map, int nt, int yday, float slat, float elat);
int rad_map_final_c(double ***rad_fract_map_final, double **rad_fract_map, int nt, int yday);
int set_min_max_hour_c(double *radfrac, double *tmin_hour, double *tmax_hour, int nx);

void set_max_min_hour(double *hourlyrad,  int tmaxhour, int tminhour); // weg?
void HourlyT_c(int nrec,
             double TmaxHour, 
             double Tmax, 
             double TminHour,
             double Tmin, 
             double *Tair);

double svp(double temp);
double sh2vp_c(double q, double t, double p);