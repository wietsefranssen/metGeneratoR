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

set_t_minmax_hour <- function(hourlyrad) {
  # /****************************************************************************
  #   set_max_min_hour
  # 
  # This function estimates the times of minimum and maximum temperature for
  # each day of the simulation, based on the hourly cycle of incoming solar
  # radiation.
  # 
  # Modifications
  # 
  # 1999-Aug-19 Modified to function in polar regions where daylight or
  # darkness may last for 24 hours.				BN
  # 2006-Oct-26 Shift tminhour and tmaxhour if necessary to remain within
  # the current day.						TJB
  # 2011-Nov-04 Changed algorithm because previous algorithm had bugs in
  # identifying the times of sunrise and sunset in some cases.
  # The new algorithm relies on the assumption that the
  # hourlyrad array is referenced to local time, i.e. hour 0 is
  # halfway between the previous day's sunset and the current
  # 	      day's sunrise.  Another assumption is that the hourlyrad
  # array begins on hour 0 of the first day.			TJB
  # 2012-Jan-28 Added logic to prevent overstepping bounds of hourlyrad
  # array.							TJB
  # ****************************************************************************/
  
  risehour = -999;
  sethour = -999;
  for (hour in 0:11) {
    if (hour == 0) {
      hourprevhour <- 24
    } else {
      hourprevhour <- hour - 1
    }
    if (hourlyrad[(hour+1)] > 0 && (hour == 0 || hourlyrad[(hourprevhour+1)] <= 0))
      risehour = hour;
  }
  for (hour in 12:23) {
    if (hour == 1) {
      hourprevhour <- 24
    } else {
      hourprevhour <- hour - 1
    }
    if (hourlyrad[(hour+1)] <= 0 && hourlyrad[(hourprevhour+1)] > 0)
      sethour = hour;
  }
  if (sethour == -999) sethour = 23;
  if (risehour >= 0 && sethour >= 0) {
    tmaxhour = 0.67 * (sethour - risehour) + risehour;
    tminhour = risehour - 1;
  } else {
    # /* arbitrarily set the min and max times to 2am and 2pm */
    tminhour = 2;
    tmaxhour = 14;
  }
  
  return(c(tminhour,tmaxhour))
}


set_max_min_hour <- function(hourlyrad) {
  # /****************************************************************************
  #   set_max_min_hour
  # 
  # This function estimates the times of minimum and maximum temperature for
  # each day of the simulation, based on the hourly cycle of incoming solar
  # radiation.
  # 
  # Modifications
  # 
  # 1999-Aug-19 Modified to function in polar regions where daylight or
  # darkness may last for 24 hours.				BN
  # 2006-Oct-26 Shift tminhour and tmaxhour if necessary to remain within
  # the current day.						TJB
  # 2011-Nov-04 Changed algorithm because previous algorithm had bugs in
  # identifying the times of sunrise and sunset in some cases.
  # The new algorithm relies on the assumption that the
  # hourlyrad array is referenced to local time, i.e. hour 0 is
  # halfway between the previous day's sunset and the current
  # 	      day's sunrise.  Another assumption is that the hourlyrad
  # array begins on hour 0 of the first day.			TJB
  # 2012-Jan-28 Added logic to prevent overstepping bounds of hourlyrad
  # array.							TJB
  # ****************************************************************************/
  
  # int risehour;
  # int sethour;
  # int hour;
  ndays<-366
  # hourlyrad<-hourly_rad_fract_data[,,1,1]
  tmaxhour <- array(NA,dim=24)
  tminhour <- array(NA,dim=24)
  for (i in 1:ndays) {
    # i<-1
    # hour<-1
    risehour = -999;
    sethour = -999;
    for (hour in 0:11) {
      if (hour == 0) {
        dayprevhour <- i - 1
        hourprevhour <- 24
      } else {
        dayprevhour <- i
        hourprevhour <- hour - 1
      }
      if (dayprevhour < 1) dayprevhour <- 366
      if (hourlyrad[i,(hour+1)] > 0 && ((i-1)*24+hour == 0 || hourlyrad[dayprevhour,(hourprevhour+1)] <= 0))
        risehour = hour;
    }
    for (hour in 12:23) {
      if (hour == 1) {
        dayprevhour <- i - 1
        hourprevhour <- 24
      } else {
        dayprevhour <- i
        hourprevhour <- hour - 1
      }
      if (dayprevhour < 1) dayprevhour <- 366
      if (hourlyrad[i,(hour+1)] <= 0 && hourlyrad[dayprevhour,(hourprevhour+1)] > 0)
        sethour = hour;
    }
    if (i == ndays && sethour == -999) sethour = 23;
    if (risehour >= 0 && sethour >= 0) {
      tmaxhour[i] = 0.67 * (sethour - risehour) + risehour;
      tminhour[i] = risehour - 1;
    }
    else {
      # /* arbitrarily set the min and max times to 2am and 2pm */
      tminhour[i] = 2;
      tmaxhour[i] = 14;
    }
    
  }
  max_min_hour<-NULL
  max_min_hour$tminhour<-tminhour
  max_min_hour$tmaxhour<-tmaxhour
  return(max_min_hour)
  
}
