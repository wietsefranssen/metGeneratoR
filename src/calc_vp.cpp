#include <math.h>
// below is based on the 'humidity' package

// "Clausius-Clapeyron"
//// [[Rcpp::export]]
float SVP_ClaCla(float t) // input in Celsius
{
  float T0 = 273.15;
  float Es_T0 = 6.11;
  float L = 2500000;
  float Rw = 461.52;
  
  float Es;
  Es = Es_T0 * exp((L/Rw) * (1/T0 - 1/(t + T0)));
  return(Es); 
}

// "Murray"
//// [[Rcpp::export]]
float SVP_Murray(float t) // input in Celsius
{
  float T0 = 273.15;
  float a, b;

  if (t < 0) {
    a = 21.8745584;
    b = 7.66;
  }
  else {
    a = 17.2693882;
    b = 35.86;
  }
  
    
  float Es;
  Es = 6.1078 * exp(a * t/((t + T0) - b));
  return(Es); 
}


//// [[Rcpp::export]]
float SH2RH(float q, float t, float p)
{
  float e, Es, psi;
  
  e = q * p/(0.622 + 0.378 * q);
  Es = SVP_ClaCla(t);
  psi = e/Es;
  
  return(psi);
}


float WVP2(float psi, float Es)
{
  float e = psi * Es;
  
  return(e);
}

//// [[Rcpp::export]]
float WVP(float psi, float t)
{
  float Es = SVP_ClaCla(t);
  float e = psi * Es;
  
  return(e);
}

//// [[Rcpp::export]]
float sh2vp_c(float q, float t, float p)
{
  float rh = SH2RH(q, t, p);
  float vp = WVP(rh, t);
  
  return(vp);
}
