HourlyT <- function(TmaxHour, Tmax,
                    TminHour, Tmin) {
  # double *x;
  # double *Tyc1;
  # double *yc2;
  # double *yc3;
  # double *yc4;
  # int i;
  # int j;
  # int n;
  # int hour;
  # int nsteps;
  HOURSPERDAY<-24
  nsteps = 24;
  
  
  
  # TmaxHour<-22
  # Tmax<-30
  # TminHour<-8
  # Tmin<-4
  # n     = ndays*2+2;
  n<- 4
  x<-Tyc1<-yc2<-yc3<-yc4<-array(0,dim=4)
  # /* First fill the x vector with the times for Tmin and Tmax, and fill the 
  # Tyc1 with the corresponding temperature and humidity values */
  # for (i = 0, j = 1, hour = 0; i < ndays; i++, hour += HOURSPERDAY) {
  # for (i = 0, j = 1, hour = 0; i < ndays; i++, hour += HOURSPERDAY) {
  # j<-1
  hour<-0
  if (TminHour < TmaxHour) {
    x[2]       = TminHour + hour
    Tyc1[2]  = Tmin
    x[3]       = TmaxHour + hour
    Tyc1[3]  = Tmax
  } else {
    x[2]       = TmaxHour + hour
    Tyc1[2]  = Tmax
    x[3]       = TminHour + hour
    Tyc1[3]  = Tmin
  } 
  
  
  # /* To "tie" down the first and last values, repeat those */
  x[1] <- x[3] - HOURSPERDAY;
  Tyc1[1] <- Tyc1[3];
  x[4] <- x[2] + HOURSPERDAY;
  # Tyc1[2] <- Tyc1[4];
  Tyc1[4] <- Tyc1[2];
  
  # /* we want to preserve maxima and minima, so we require that the first 
  # derivative at these points is zero */
  for (i in 1:4) {
    yc2[i] = 0.;
  }
  
  # /* calculate the coefficients for the splines for the temperature */
  # hermite(n, x, Tyc1, yc2, yc3, yc4);
  for (i in 1:n-1) {
    # for (i = 0; i < n-1; i++) {
    dx = x[i+1] - x[i];
    divdf1 = (Tyc1[i+1] - Tyc1[i])/dx;
    # divdf1 = (yc1[i+1] - yc1[i])/dx;
    divdf3 = yc2[i] + yc2[i+1] - 2 * divdf1;
    yc3[i] = (divdf1 - yc2[i] - divdf3)/dx;
    yc4[i] = divdf3/(dx*dx);
  }
  
  # /* interpolate the temperatures */
  # for (i = 0, hour = 0; i < nsteps; i++, hour += Dt) {
  hour<-0
  Tair<-array(0,dim=24)
  for (i in 1:nsteps) {
    # i<-2
    # Tair[i] = hermint(hour, n, x, Tyc1, yc2, yc3, yc4);
    # hermint(xbar, n, double *x, double *yc1, double *yc2, 
    #         double *yc3, double *yc4)
    
    klo=0;
    khi=n-1;
    while (khi-klo > 1) {
      k=bitwShiftR((khi+klo), 1)
      # k=(khi+klo) >> 1;
      if (x[(k+1)] > hour) {
        khi=k 
      } else {
        klo=k
      }
    }
    
    dx = hour - x[klo+1];
    Tair[i] = Tyc1[klo+1] + dx * (yc2[klo+1] + dx * (yc3[klo+1] + dx * yc4[klo+1]));
    # return result;
    hour<-hour+1
  }
  # print(Tair)
  # 
  # hour <- i-1
  return(Tair)
}

# /**************************************************************************/

hermite <- function(n, x, yc1, yc2, yc3, yc4) {
  # int i;
  # double dx;
  # double divdf1;
  # double divdf3;
  
  for (i in 1:n-1) {
    # for (i = 0; i < n-1; i++) {
    dx = x[i+1] - x[i];
    divdf1 = (yc1[i+1] - yc1[i])/dx;
    divdf3 = yc2[i] + yc2[i+1] - 2 * divdf1;
    yc3[i] = (divdf1 - yc2[i] - divdf3)/dx;
    yc4[i] = divdf3/(dx*dx);
  }
}
