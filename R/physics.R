# require(compiler)
# 
# # source("./R/constants.R")
# 
# # cnst$SW_RAD_DT<-30
# #cnst$SW_RAD_DT<-30*10
# 
# solar_geom_old <- function(elev = elev, lat = lat, lr = lr) {
#   # """
#   # Flat earth assumption
#   #
#   # Parameters
#   # ----------
#   # elev:
#   #     Elevation in meters
#   # lat:
#   #     Latitude in decimal format
#   # lr:
#   #     Lapse rate in K/m
#   #
#   # Returns
#   # -------
#   # sg:
#   #     (tiny_rad_fract, daylength, flat_potrad, tt_max0)
#   # """
# 
#   # optical airmass by degrees
#   OPTAM <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
#              4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
#              8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
#   dayperyear <- ceiling(cnst$DAYS_PER_YEAR)
#   tt_max0 <- rep(0, dayperyear)
#   daylength <-  rep(0, dayperyear)
#   flat_potrad <-  rep(0, dayperyear)
# 
#   # Calculate pressure ratio as a function of elevation
#   t1 <- 1.0 - (lr * elev) / cnst$T_STD
#   t2 <- cnst$G_STD / (lr * (cnst$R / cnst$MA))
#   trans <- cnst$TBASE^(t1^t2)
# 
#   # Translate lat to rad
#   lat <- min(max(lat * cnst$RAD_PER_DEG, -pi / 2.), pi / 2.0)
#   coslat <- cos(lat)
#   sinlat <- sin(lat)
# 
#   # Sub-daily time step and angular step
#   dt <- cnst$SW_RAD_DT
#   dh <- dt / cnst$SEC_PER_RAD
# 
#   # Allocate the radiation arrays
#   tiny_step_per_day <- cnst$SEC_PER_DAY / cnst$SW_RAD_DT
#   tiny_rad_fract <- array(0, dim = c(dayperyear, tiny_step_per_day))
#   for (i in 1:dayperyear) {
#     # Declination and quantities of interest
#     decl <- cnst$MIN_DECL * cos(((i-1) + cnst$DAYS_OFF) * cnst$RAD_PER_DAY)
#     cosdecl <- cos(decl)
#     sindecl <- sin(decl)
# 
#     # calculate daylength as a function of lat and decl
#     cosegeom <- coslat * cosdecl
#     sinegeom <- sinlat * sindecl
#     coshss <- min(max(-sinegeom / cosegeom, -1), 1)
#     hss <- acos(coshss)
#     daylength[i] <- min(2.0 * hss * cnst$SEC_PER_RAD, cnst$SEC_PER_DAY)
#     # Extraterrestrial radiation perpendicular to beam,
#     # total over the timestep (J)
#     dir_beam_topa <- (1368.0 + 45.5 * sin((2.0 * pi * (i-1) / cnst$DAYS_PER_YEAR) + 1.7)) * dt
#     sum_trans <- 0
#     sum_flat_potrad <- 0
#     # Set up angular calculations
#     for (h in seq(from = -hss, to = hss, by = dh)) {
#       # Cosine of the hour angle and solar zenith angle
#       cosh <- cos(h)
#       cza <- cosegeom * cosh + sinegeom
#       if (cza > 0) {
#         # When sun is above flat horizon do flat-surface
#         # calculations to determine daily total transmittance
#         # and save potential radiation for calculation of
#         # diffuse portion
#         dir_flat_topa <- dir_beam_topa * cza
#         am <- 1.0 / (cza + 0.0000001)
#         if (am > 2.9) {
#           ami <- min(max((acos(cza) / cnst$RAD_PER_DEG) - 69, 0), 20)
#           # print(ami)
#           am <- OPTAM[ami+1]
#           # print(am)
#         }
#         sum_trans <- sum_trans + ((trans^am) * dir_flat_topa)
#         sum_flat_potrad <- sum_flat_potrad + dir_flat_topa
#       } else {
#         # Sun not above horizon
#         dir_flat_topa <- 0
#       }
#       tinystep <- min(max((12 * cnst$SEC_PER_HOUR + h * cnst$SEC_PER_RAD) / dt, 0),tiny_step_per_day - 1)
#       tiny_rad_fract[i,tinystep] <- dir_flat_topa
#     }
#     if (daylength[i] && sum_flat_potrad > 0) {
#       tiny_rad_fract[i,] <- tiny_rad_fract[i,] / sum_flat_potrad
#     }
#     if (daylength[i]) {
#       # Transmittance and potential radiation
#       # averaged over daylength
#       # print(sum_trans )
#       # print( sum_flat_potrad)
#       # print(sum_trans / sum_flat_potrad)
#       tt_max0[i] <- sum_trans / sum_flat_potrad
#       flat_potrad[i] <- sum_flat_potrad / daylength[i]
#     } else {
#       # No daytime - no radiation
#       tt_max0[i] <- 0.
#       flat_potrad[i] <- 0.
#     }
#   }
#   tt_max0[dayperyear] <- tt_max0[dayperyear - 1]
#   flat_potrad[dayperyear] <- flat_potrad[dayperyear - 1]
#   daylength[dayperyear] <- daylength[dayperyear - 1]
#   tiny_rad_fract[dayperyear] <- tiny_rad_fract[dayperyear - 1]
#   # return tiny_rad_fract, daylength, flat_potrad, tt_max0
#   all<-NULL
#   all$tt_max0<-tt_max0
#   all$flat_potrad<-flat_potrad
#   all$daylength<-daylength
#   all$tiny_rad_fract<-tiny_rad_fract
#   return(all)
# }
# 
# solar_geom <- function(lat = lat) {
#   # """
#   # Flat earth assumption
#   #
#   # Parameters
#   # ----------
#   # lat:
#   #     Latitude in decimal format
#   #
#   # Returns
#   # -------
#   # sg:
#   #     (tiny_rad_fract, daylength, flat_potrad, tt_max0)
#   # """
# 
#   # optical airmass by degrees
#   OPTAM <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
#              4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
#              8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
#   dayperyear <- ceiling(cnst$DAYS_PER_YEAR)
#   tt_max0 <- rep(0, dayperyear)
#   daylength <-  rep(0, dayperyear)
#   flat_potrad <-  rep(0, dayperyear)
# 
#   # Translate lat to rad
#   lat <- min(max(lat * cnst$RAD_PER_DEG, -pi / 2.), pi / 2.0)
#   coslat <- cos(lat)
#   sinlat <- sin(lat)
# 
#   # Sub-daily time step and angular step
#   dt <- cnst$SW_RAD_DT
#   dh <- dt / cnst$SEC_PER_RAD
# 
#   # Allocate the radiation arrays
#   tiny_step_per_day <- cnst$SEC_PER_DAY / cnst$SW_RAD_DT
#   tiny_rad_fract <- array(0, dim = c(dayperyear, tiny_step_per_day))
#   for (i in 1:dayperyear) {
#     # Declination and quantities of interest
#     decl <- cnst$MIN_DECL * cos(((i-1) + cnst$DAYS_OFF) * cnst$RAD_PER_DAY)
#     cosdecl <- cos(decl)
#     sindecl <- sin(decl)
# 
#     # calculate daylength as a function of lat and decl
#     cosegeom <- coslat * cosdecl
#     sinegeom <- sinlat * sindecl
#     coshss <- min(max(-sinegeom / cosegeom, -1), 1)
#     hss <- acos(coshss)
#     daylength[i] <- min(2.0 * hss * cnst$SEC_PER_RAD, cnst$SEC_PER_DAY)
#     # Extraterrestrial radiation perpendicular to beam,
#     # total over the timestep (J)
#     dir_beam_topa <- (1368.0 + 45.5 * sin((2.0 * pi * (i-1) / cnst$DAYS_PER_YEAR) + 1.7)) * dt
#     sum_trans <- 0
#     sum_flat_potrad <- 0
#     # Set up angular calculations
#     for (h in seq(from = -hss, to = hss, by = dh)) {
#       # Cosine of the hour angle and solar zenith angle
#       cosh <- cos(h)
#       cza <- cosegeom * cosh + sinegeom
#       if (cza > 0) {
#         # When sun is above flat horizon do flat-surface
#         # calculations to determine daily total transmittance
#         # and save potential radiation for calculation of
#         # diffuse portion
#         dir_flat_topa <- dir_beam_topa * cza
#         am <- 1.0 / (cza + 0.0000001)
#         if (am > 2.9) {
#           ami <- min(max((acos(cza) / cnst$RAD_PER_DEG) - 69, 0), 20)
#           # print(ami)
#           am <- OPTAM[ami+1]
#           # print(am)
#         }
#         sum_flat_potrad <- sum_flat_potrad + dir_flat_topa
#       } else {
#         # Sun not above horizon
#         dir_flat_topa <- 0
#       }
#       tinystep <- min(max((12 * cnst$SEC_PER_HOUR + h * cnst$SEC_PER_RAD) / dt, 0), tiny_step_per_day - 1)
#       tiny_rad_fract[i,tinystep] <- dir_flat_topa
#     }
#     if (daylength[i] && sum_flat_potrad > 0) {
#       tiny_rad_fract[i,] <- tiny_rad_fract[i,] / sum_flat_potrad
#     }
#     if (daylength[i]) {
#       # Transmittance and potential radiation
#       # averaged over daylength
#       flat_potrad[i] <- sum_flat_potrad / daylength[i]
#     } else {
#       # No daytime - no radiation
#       tt_max0[i] <- 0.
#       flat_potrad[i] <- 0.
#     }
#   }
#   flat_potrad[dayperyear] <- flat_potrad[dayperyear - 1]
#   daylength[dayperyear] <- daylength[dayperyear - 1]
#   tiny_rad_fract[dayperyear] <- tiny_rad_fract[dayperyear - 1]
# 
#   all<-NULL
#   all$flat_potrad<-flat_potrad
#   all$daylength<-daylength
#   all$tiny_rad_fract<-tiny_rad_fract
#   return(all)
# }
# 
# solar_geom_tt_max0 <- function(elev = elev, lat = lat, lr = lr) {
#   # """
#   # Flat earth assumption
#   #
#   # Parameters
#   # ----------
#   # elev:
#   #     Elevation in meters
#   # lat:
#   #     Latitude in decimal format
#   # lr:
#   #     Lapse rate in K/m
#   #
#   # Returns
#   # -------
#   # sg:
#   #     (tiny_rad_fract, daylength, flat_potrad, tt_max0)
#   # """
# 
#   # optical airmass by degrees
#   OPTAM <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
#              4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
#              8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
#   dayperyear <- ceiling(cnst$DAYS_PER_YEAR)
#   tt_max0 <- rep(0, dayperyear)
#   daylength <-  rep(0, dayperyear)
#   flat_potrad <-  rep(0, dayperyear)
# 
#   # Calculate pressure ratio as a function of elevation
#   t1 <- 1.0 - (lr * elev) / cnst$T_STD
#   t2 <- cnst$G_STD / (lr * (cnst$R / cnst$MA))
#   trans <- cnst$TBASE^(t1^t2)
# 
#   # Translate lat to rad
#   lat <- min(max(lat * cnst$RAD_PER_DEG, -pi / 2.), pi / 2.0)
#   coslat <- cos(lat)
#   sinlat <- sin(lat)
# 
#   # Sub-daily time step and angular step
#   dt <- cnst$SW_RAD_DT
#   dh <- dt / cnst$SEC_PER_RAD
# 
#   # Allocate the radiation arrays
#   for (i in 1:dayperyear) {
#     # Declination and quantities of interest
#     decl <- cnst$MIN_DECL * cos(((i-1) + cnst$DAYS_OFF) * cnst$RAD_PER_DAY)
#     cosdecl <- cos(decl)
#     sindecl <- sin(decl)
# 
#     # calculate daylength as a function of lat and decl
#     cosegeom <- coslat * cosdecl
#     sinegeom <- sinlat * sindecl
#     coshss <- min(max(-sinegeom / cosegeom, -1), 1)
#     hss <- acos(coshss)
#     daylength[i] <- min(2.0 * hss * cnst$SEC_PER_RAD, cnst$SEC_PER_DAY)
#     # Extraterrestrial radiation perpendicular to beam,
#     # total over the timestep (J)
#     dir_beam_topa <- (1368.0 + 45.5 * sin((2.0 * pi * (i-1) / cnst$DAYS_PER_YEAR) + 1.7)) * dt
#     sum_trans <- 0
#     sum_flat_potrad <- 0
#     # Set up angular calculations
#     for (h in seq(from = -hss, to = hss, by = dh)) {
#       # Cosine of the hour angle and solar zenith angle
#       cosh <- cos(h)
#       cza <- cosegeom * cosh + sinegeom
#       if (cza > 0) {
#         # When sun is above flat horizon do flat-surface
#         # calculations to determine daily total transmittance
#         # and save potential radiation for calculation of
#         # diffuse portion
#         dir_flat_topa <- dir_beam_topa * cza
#         am <- 1.0 / (cza + 0.0000001)
#         if (am > 2.9) {
#           ami <- min(max((acos(cza) / cnst$RAD_PER_DEG) - 69, 0), 20)
#           # print(ami)
#           am <- OPTAM[ami+1]
#           # print(am)
#         }
#         sum_trans <- sum_trans + ((trans^am) * dir_flat_topa)
#         sum_flat_potrad <- sum_flat_potrad + dir_flat_topa
#       } else {
#         # Sun not above horizon
#         dir_flat_topa <- 0
#       }
#     }
#     if (daylength[i]) {
#       # Transmittance and potential radiation
#       # averaged over daylength
#       tt_max0[i] <- sum_trans / sum_flat_potrad
#     } else {
#       # No daytime - no radiation
#       tt_max0[i] <- 0.
#     }
#   }
#   tt_max0[dayperyear] <- tt_max0[dayperyear - 1]
#   return(tt_max0)
# }
# 
# csolar_geom_old <- cmpfun(solar_geom_old)
# csolar_geom <- cmpfun(solar_geom)
# csolar_geom_tt_max0 <- cmpfun(solar_geom_tt_max0)
# 
