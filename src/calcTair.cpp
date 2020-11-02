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
             float *x, 
             float *yc1, 
             float *yc2, 
             float *yc3, 
             float *yc4)
{
  int i;
  float dx;
  float divdf1;
  float divdf3;
  
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
float hermint(float xbar, int n, float *x, float *yc1, float *yc2, 
               float *yc3, float *yc4)
{
  int klo,khi,k;
  float dx;
  float result;
  
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
               float TminHour,
               float Tmin, 
               float TmaxHour, 
               float Tmax, 
               float *Tair) 
{
  float *x;
  float *Tyc1;
  float *yc2;
  float *yc3;
  float *yc4;
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
  x     = (float *) calloc(n, sizeof(float));
  Tyc1  = (float *) calloc(n, sizeof(float));
  yc2   = (float *) calloc(n, sizeof(float));
  yc3   = (float *) calloc(n, sizeof(float));
  yc4   = (float *) calloc(n, sizeof(float));

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

