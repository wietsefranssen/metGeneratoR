#include <math.h>
// below is based on the 'humidity' package

// "Clausius-Clapeyron"
//// [[Rcpp::export]]
double SVP_ClaCla(double t) // input in Celsius
{
  double T0 = 273.15;
  double Es_T0 = 6.11;
  double L = 2500000;
  double Rw = 461.52;
  
  double Es;
  Es = Es_T0 * exp((L/Rw) * (1/T0 - 1/(t + T0)));
  return(Es); 
}

// "Murray"
//// [[Rcpp::export]]
double SVP_Murray(double t) // input in Celsius
{
  double T0 = 273.15;
  double a, b;

  if (t < 0) {
    a = 21.8745584;
    b = 7.66;
  }
  else {
    a = 17.2693882;
    b = 35.86;
  }
  
    
  double Es;
  Es = 6.1078 * exp(a * t/((t + T0) - b));
  return(Es); 
}


//// [[Rcpp::export]]
double SH2RH(double q, double t, double p)
{
  double e, Es, psi;
  
  e = q * p/(0.622 + 0.378 * q);
  Es = SVP_ClaCla(t);
  psi = e/Es;
  
  return(psi);
}


double WVP2(double psi, double Es)
{
  double e = psi * Es;
  
  return(e);
}

//// [[Rcpp::export]]
double WVP(double psi, double t)
{
  double Es = SVP_ClaCla(t);
  double e = psi * Es;
  
  return(e);
}

//// [[Rcpp::export]]
double sh2vp_c(double q, double t, double p)
{
  double rh = SH2RH(q, t, p);
  double vp = WVP(rh, t);
  
  return(vp);
}
