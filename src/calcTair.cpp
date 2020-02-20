#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/****************************************************************************
 Subroutines developed by Bart Nijssen to estimate the daily temperature
 cycle from maximum and minimum daily temperature measurements.  
 
 Modifications:
 June 23, 1998 by Keith Cherkauer to be run within the VIC-NL model.
 *****************************************************************************/

/****************************************************************************/
/*				    hermite                                 */
/****************************************************************************/
/* calculate the coefficients for the Hermite polynomials */
void hermite(int n, 
             double *x, 
             double *yc1, 
             double *yc2, 
             double *yc3, 
             double *yc4)
{
  int i;
  double dx;
  double divdf1;
  double divdf3;
  
  for (i = 0; i < n-1; i++) {
    dx = x[i+1] - x[i];
    divdf1 = (yc1[i+1] - yc1[i])/dx;
    divdf3 = yc2[i] + yc2[i+1] - 2 * divdf1;
    yc3[i] = (divdf1 - yc2[i] - divdf3)/dx;
    yc4[i] = divdf3/(dx*dx);
  }
}

/**************************************************************************/
/*				    hermint                               */
/**************************************************************************/
/* use the Hermite polynomials, to find the interpolation function value at 
 xbar */
double hermint(double xbar, int n, double *x, double *yc1, double *yc2, 
               double *yc3, double *yc4)
{
  int klo,khi,k;
  double dx;
  double result;
  
  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (x[k] > xbar) khi=k;
    else klo=k;
  }
  
  dx = xbar - x[klo];
  result = yc1[klo] + dx * (yc2[klo] + dx * (yc3[klo] + dx * yc4[klo]));
  return result;
}

/****************************************************************************/
/*				    HourlyT                                 */
/****************************************************************************/
void HourlyT_c(int nrec,
               double TminHour,
               double Tmin, 
               double TmaxHour, 
               double Tmax, 
               double *Tair) 
{
  double *x;
  double *Tyc1;
  double *yc2;
  double *yc3;
  double *yc4;
  int i;
  int n;
  int hour;
  int nsteps;
  
  int HOURSPERDAY = nrec;
  TminHour = TminHour / (24/nrec);
  TmaxHour = TmaxHour / (24/nrec);
  // int HOURSPERDAY = 24;
  // nsteps = HOURSPERDAY/Dt * ndays;
  nsteps = nrec;
  
  n     = 2+2;
  x     = (double *) calloc(n, sizeof(double));
  Tyc1  = (double *) calloc(n, sizeof(double));
  yc2   = (double *) calloc(n, sizeof(double));
  yc3   = (double *) calloc(n, sizeof(double));
  yc4   = (double *) calloc(n, sizeof(double));
  
  /* First fill the x vector with the times for Tmin and Tmax, and fill the 
   Tyc1 with the corresponding temperature and humidity values */
  hour = 0;
  if (TminHour < TmaxHour) {
    x[1]       = TminHour + hour;
    Tyc1[1]  = Tmin;
    x[2]       = TmaxHour + hour;
    Tyc1[2]  = Tmax;
  }
  else {
    // printf("TminHour > TmaxHour");
    x[1]       = TmaxHour + hour;
    Tyc1[1]  = Tmax;
    x[2]       = TminHour + hour;
    Tyc1[2]  = Tmin;
  } 
  
  /* To "tie" down the first and last values, repeat those */
  x[0] = x[2] - HOURSPERDAY;
  Tyc1[0] = Tyc1[2];
  x[3] = x[1] + HOURSPERDAY;
  Tyc1[3] = Tyc1[1];
  
  /* we want to preserve maxima and minima, so we require that the first 
   derivative at these points is zero */
  for (i = 0; i < n; i++)
    yc2[i] = 0.;
  
  /* calculate the coefficients for the splines for the temperature */
  // hermite(n, x, Tyc1, yc2, yc3, yc4);
  // for (i in 1:n-1) {
  float dx, divdf1,divdf3;
  for (i = 0; i < (n-1); i++)
  {
    dx = x[i+1] - x[i];
    divdf1 = (Tyc1[i+1] - Tyc1[i])/dx;
    divdf3 = yc2[i] + yc2[i+1] - 2 * divdf1;
    yc3[i] = (divdf1 - yc2[i] - divdf3)/dx;
    yc4[i] = divdf3/(dx*dx);
  }
  
  /* interpolate the temperatures */
  // for (hour = 0; hour < nsteps; hour++) {
  //   Tair[i] = hermint(hour, n, x, Tyc1, yc2, yc3, yc4);
  // }
  // for (hour = 0; hour < nsteps; hour++) {
  //   Tair[hour] = 0;
  // }
  
  int klo,khi,k;
  hour = 0;
  for (i = 0; i < nsteps; i++) {
    
    klo=0;
    khi=n-1;
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      
      if (x[k] > hour) {
        khi=k;
      } else {
        klo=k;
      }
    }
    
    dx = hour - x[klo];
    Tair[i] = Tyc1[klo] + dx * (yc2[klo] + dx * (yc3[klo] + dx * yc4[klo]));
    hour++;
  }
  
  free(x);   
  free(Tyc1);
  free(yc2);
  free(yc3);
  free(yc4);
  
  return;
}

void set_max_min_hour(double *hourlyrad, 
                      int tmaxhour,
                      int tminhour)
  /****************************************************************************
   set_max_min_hour
   
   This function estimates the times of minimum and maximum temperature for
   each day of the simulation, based on the hourly cycle of incoming solar
   radiation.
   
   Modifications
   
   1999-Aug-19 Modified to function in polar regions where daylight or
   darkness may last for 24 hours.				BN
   2006-Oct-26 Shift tminhour and tmaxhour if necessary to remain within
   the current day.						TJB
   2011-Nov-04 Changed algorithm because previous algorithm had bugs in
   identifying the times of sunrise and sunset in some cases.
   The new algorithm relies on the assumption that the
   hourlyrad array is referenced to local time, i.e. hour 0 is
   halfway between the previous day's sunset and the current
   day's sunrise.  Another assumption is that the hourlyrad
   array begins on hour 0 of the first day.			TJB
   2012-Jan-28 Added logic to prevent overstepping bounds of hourlyrad
   array.							TJB
   ****************************************************************************/
{
  int risehour;
  int sethour;
  int hour;
  
  risehour = -999;
  sethour = -999;
  for (hour = 0; hour < 12; hour++) {
    printf("rad: %f\n", hourlyrad[hour]);
    if (hourlyrad[hour] > 0 && (hour == 0 || hourlyrad[hour-1] <= 0))
      risehour = hour;
  }
  for (hour = 12; hour < 24; hour++) {
    printf("rad: %f\n", hourlyrad[hour]);
    if (hourlyrad[hour] <= 0 && hourlyrad[hour-1] > 0)
      sethour = hour;
  }
  if (sethour == -999) sethour = 23;
  if (risehour >= 0 && sethour >= 0) {
    tmaxhour = 0.67 * (sethour - risehour) + risehour;
    tminhour = risehour - 1;
  }
  else {
    /* arbitrarily set the min and max times to 2am and 2pm */
    tminhour = 2;
    tmaxhour = 14;
  }
  
  printf("min: %d, max: %d\n", tminhour, tmaxhour);
  
  
}
