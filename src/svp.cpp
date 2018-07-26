#include <math.h>

double svp(double temp)
/**********************************************************************
  This routine computes the saturated vapor pressure using Handbook
of Hydrology eqn 4.2.2

Pressure in Pa

**********************************************************************/
{
  double SVP;
  
double A_SVP=0.61078;
double B_SVP =17.269;
double C_SVP =237.3;
  
  SVP = A_SVP * exp((B_SVP * temp)/(C_SVP+temp));
  
  if(temp<0) SVP *= 1.0 + .00972 * temp + .000042 * temp * temp;
  
  return (SVP*1000.);
}