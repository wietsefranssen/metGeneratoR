mtclim_to_vic <- function(hour_offset, yday, tiny_radfract, ctrl, 
                          mtclim_data, tskc, vp, fdir)
  # void mtclim_to_vic(double hour_offset, 
  #                    int Ndays, dmy_struct *dmy, 
  #                    double **tiny_radfract, control_struct *ctrl, 
  #                    data_struct *mtclim_data, double *tskc, double *vp, 
  #                    double *hourlyrad, double *fdir)
  # /******************************************************************************
  #   mtclim_to_vic: Store MTCLIM variables in VIC arrays.
  # 
  # Modifications:
  #   2012-Feb-16 Removed check on mtclim_data->insw for storing tinyradfract data
  # in hourlyrad array.						TJB
# ******************************************************************************/
{
  tinystepsphour <- 3600/metGen$constants$SRADDT;
  
  tiny_offset <- tinystepsphour * hour_offset
  # // s_srad = avg SW flux (W/m2) over daylight hours
  # // s_dayl = number of seconds of daylight in current day
  # // total_daily_sw = s_srad*s_dayl  (J/m2)
  # // tiny_radfrac = fraction of total daily sw falling in each SRADDT interval
  # // hourlyrad = SW flux (W/m2) over each hour = total_daily_sw * sum_over_hour(tiny_radfract) / 3600
  # //                                           = tmp_rad * sum_over_hour(tiny_radfract)
  
  # //if radiation read from input file, assume it's a 24 hours average, 
  # //else (i.e., MTCLIM calculated), assume it's a daylight period average
  hourlyrad <- array(NA, dim = (24))
  
  if (ctrl$insw) {
    tmp_rad <- mtclim_data$s_srad * 24.;
  } else {
    tmp_rad = mtclim_data$s_srad * mtclim_data$s_dayl / 3600.;
  }    
  for (j in 1:24) {
    hourlyrad[j] = 0;
    for (k in 1:tinystepsphour) {
      tinystep = ((j-1)*tinystepsphour + (k-1) - tiny_offset) + 1
      if (tinystep < 1) {
        tinystep <- tinystep + 24*tinystepsphour; 
        # print(paste0("<1 ", tinystep))
      }
      if (tinystep > 24*tinystepsphour) {
        tinystep <- tinystep - 24*tinystepsphour; 
      }
      hourlyrad[j] <- hourlyrad[j] + tiny_radfract[yday,tinystep];
    }
    hourlyrad[j] <- hourlyrad[j] * tmp_rad;
  }

  fdir = mtclim_data$s_fdir;
  tskc = mtclim_data$s_tskc;
  vp = mtclim_data$s_hum;
  return(hourlyrad)
}

