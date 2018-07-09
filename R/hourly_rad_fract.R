hourly_rad_fract2 <- function(mask, theta_l, tiny_rad_fract)
  # void mtclim_to_vic(double hour_offset, 
  #                    int Ndays, dmy_struct *dmy, 
  #                    double **tiny_rad_fract, control_struct *ctrl, 
  #                    data_struct *mtclim_data, double *tskc, double *vp, 
  #                    double *hourly_rad_fract, double *fdir)
  # /******************************************************************************
  #   mtclim_to_vic: Store MTCLIM variables in VIC arrays.
  # 
  # Modifications:
  #   2012-Feb-16 Removed check on mtclim_data->insw for storing tinyradfract data
  # in hourly_rad_fract array.						TJB
# ******************************************************************************/
{
  hourly_rad_fract <- array(0, dim = c(366, 24, length(mask$xyCoords$y)))
  
 for (day in 1:366) {
    for (k in 1:(86400/30)) {
      hour<-ceiling(k/tinystepsphour)
      
      hourly_rad_fract[day,hour,] <- hourly_rad_fract[day,hour,] + tiny_rad_fract[day, k,]
    }
  }
  
  return(hourly_rad_fract)
}


hourly_rad_fract <- function(mask, theta_l, tiny_rad_fract)
  # void mtclim_to_vic(double hour_offset, 
  #                    int Ndays, dmy_struct *dmy, 
  #                    double **tiny_rad_fract, control_struct *ctrl, 
  #                    data_struct *mtclim_data, double *tskc, double *vp, 
  #                    double *hourly_rad_fract, double *fdir)
  # /******************************************************************************
  #   mtclim_to_vic: Store MTCLIM variables in VIC arrays.
  # 
  # Modifications:
  #   2012-Feb-16 Removed check on mtclim_data->insw for storing tinyradfract data
  # in hourly_rad_fract array.						TJB
# ******************************************************************************/
{
  # tiny_rad_fract <- solar_geom$tiny_rad_fract
  
  ## Define julian days for the rest of the year based on first day of run
  yday<-array(NA, 366)
  yday_start<-metGen$derived$inYDays[1]
  yday_tmp <- yday_start
  for (i in 1:366) {
    if (yday_tmp > 366) yday_tmp <- 1
    yday[i] <- yday_tmp
    yday_tmp <- yday_tmp + 1
  }
  
  hourly_rad_fract <- array(0, dim = c(366, 24, length(mask$xyCoords$x), length(mask$xyCoords$y)))
  for (ilon in 1:length(mask$xyCoords$x)) {
    ## Calculate offset longitude
    hour_offset = (theta_l-mask$xyCoords$x[ilon])*24/360;
    if (hour_offset < 0) {
      hour_offset_int = floor(hour_offset-0.5);
    } else {
      hour_offset_int = floor(hour_offset+0.5);
    }
    hour_offset <- hour_offset - hour_offset_int #// hour_offset is now the distance from the center of local time zone
    
    tinystepsphour <- 3600/metGen$constants$SRADDT;
    
    tiny_offset <- tinystepsphour * hour_offset
    # // s_srad = avg SW flux (W/m2) over daylight hours
    # // s_dayl = number of seconds of daylight in current day
    # // total_daily_sw = s_srad*s_dayl  (J/m2)
    # // tiny_radfrac = fraction of total daily sw falling in each SRADDT interval
    # // hourly_rad_fract = SW flux (W/m2) over each hour = total_daily_sw * sum_over_hour(tiny_rad_fract) / 3600
    # //                                           = tmp_rad * sum_over_hour(tiny_rad_fract)
    
    # //if radiation read from input file, assume it's a 24 hours average, 
    # //else (i.e., MTCLIM calculated), assume it's a daylight period average
    
    
    for (day in 1:366) {
      for (j in 1:24) {
        hourly_rad_fract[j] = 0;
        for (k in 1:tinystepsphour) {
          tinystep = ((j-1)*tinystepsphour + (k-1) - tiny_offset) + 1
          # print(tinystep)
          if (tinystep < 1) {
            tinystep <- tinystep + 24*tinystepsphour; 
          }
          if (tinystep > 24*tinystepsphour) {
            tinystep <- tinystep - 24*tinystepsphour; 
          }
          hourly_rad_fract[day,j,ilon,] <- hourly_rad_fract[day,j,ilon,] + tiny_rad_fract[yday[day],tinystep,]
        }
      }
    }
  }
  return(hourly_rad_fract)
}
