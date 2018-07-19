#include <math.h>
#include <Rcpp.h>
#include "Rdefines.h" // R memory io
#include "Rmath.h"    // R math functions

/*** SubRoutine Prototypes ***/
int solar_geom_c(double*, float, int , int);
int rad_fract_lats_c(double **rad_fract_map, int nt, int yday);
int rad_map_final_c(double ***rad_fract_map_final, double **rad_fract_map, int nt, int yday);

void set_max_min_hour(double *hourlyrad, 
                      int tmaxhour,
                      int tminhour);