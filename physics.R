require(compiler)
## Cleanup
rm(list=ls(all=TRUE))

setwd("/home/wietse/Documents/WORKDIRS/RProj/disaggregateR/")

elev<-500
lat<-50.25
source("./constants.R")

cnst$SW_RAD_DT<-30
cnst$SW_RAD_DT<-30*10

solar_geom <- function(elev = elev, lat = lat, lr = lr) {
  # """
  # Flat earth assumption
  #
  # Parameters
  # ----------
  # elev:
  #     Elevation in meters
  # lat:
  #     Latitude in decimal format
  # lr:
  #     Lapse rate in K/m
  #
  # Returns
  # -------
  # sg:
  #     (tiny_rad_fract, daylength, flat_potrad, tt_max0)
  # """
  
  # optical airmass by degrees
  OPTAM <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
             4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
             8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
  dayperyear <- ceiling(cnst$DAYS_PER_YEAR)
  tt_max0 <- rep(0, dayperyear)
  daylength <-  rep(0, dayperyear)
  flat_potrad <-  rep(0, dayperyear)
  
  # Calculate pressure ratio as a function of elevation
  t1 <- 1.0 - (lr * elev) / cnst$T_STD
  t2 <- cnst$G_STD / (lr * (cnst$R / cnst$MA))
  trans <- cnst$TBASE^(t1^t2)
  
  # Translate lat to rad
  lat <- min(max(lat * cnst$RAD_PER_DEG, -pi / 2.), pi / 2.0)
  coslat <- cos(lat)
  sinlat <- sin(lat)
  
  # Sub-daily time step and angular step
  dt <- cnst$SW_RAD_DT
  dh <- dt / cnst$SEC_PER_RAD
  
  # Allocate the radiation arrays
  tiny_step_per_day <- cnst$SEC_PER_DAY / cnst$SW_RAD_DT
  tiny_rad_fract <- array(0, dim = c(dayperyear, tiny_step_per_day))
  for (i in 1:dayperyear) {
    # Declination and quantities of interest
    decl <- cnst$MIN_DECL * cos(((i-1) + cnst$DAYS_OFF) * cnst$RAD_PER_DAY)
    cosdecl <- cos(decl)
    sindecl <- sin(decl)
    
    # calculate daylength as a function of lat and decl
    cosegeom <- coslat * cosdecl
    sinegeom <- sinlat * sindecl
    coshss <- min(max(-sinegeom / cosegeom, -1), 1)
    hss <- acos(coshss)
    daylength[i] <- min(2.0 * hss * cnst$SEC_PER_RAD, cnst$SEC_PER_DAY)
    # Extraterrestrial radiation perpendicular to beam,
    # total over the timestep (J)
    dir_beam_topa <- (1368.0 + 45.5 * sin((2.0 * pi * (i-1) / cnst$DAYS_PER_YEAR) + 1.7)) * dt
    sum_trans <- 0
    sum_flat_potrad <- 0
    # Set up angular calculations
    for (h in seq(from = -hss, to = hss, by = dh)) {
      # Cosine of the hour angle and solar zenith angle
      cosh <- cos(h)
      cza <- cosegeom * cosh + sinegeom
      if (cza > 0) {
        # When sun is above flat horizon do flat-surface
        # calculations to determine daily total transmittance
        # and save potential radiation for calculation of
        # diffuse portion
        dir_flat_topa <- dir_beam_topa * cza
        am <- 1.0 / (cza + 0.0000001)
        if (am > 2.9) {
          ami <- min(max((acos(cza) / cnst$RAD_PER_DEG) - 69, 0), 20)
          # print(ami)
          am <- OPTAM[ami+1]
          # print(am)
        }
        sum_trans <- sum_trans + ((trans^am) * dir_flat_topa)
        sum_flat_potrad <- sum_flat_potrad + dir_flat_topa
      } else {
        # Sun not above horizon
        dir_flat_topa <- 0
      }
      tinystep <- min(max((12 * cnst$SEC_PER_HOUR + h * cnst$SEC_PER_RAD) / dt, 0),tiny_step_per_day - 1)
      tiny_rad_fract[i,tinystep] <- dir_flat_topa
    }
    if (daylength[i] && sum_flat_potrad > 0) {
      tiny_rad_fract[i,] <- tiny_rad_fract[i,] / sum_flat_potrad
    }
    if (daylength[i]) {
      # Transmittance and potential radiation
      # averaged over daylength
      # print(sum_trans )
      # print( sum_flat_potrad)
      # print(sum_trans / sum_flat_potrad)
      tt_max0[i] <- sum_trans / sum_flat_potrad
      flat_potrad[i] <- sum_flat_potrad / daylength[i]
    } else {
      # No daytime - no radiation
      tt_max0[i] <- 0.
      flat_potrad[i] <- 0.
    }
  }
  tt_max0[dayperyear] <- tt_max0[dayperyear - 1]
  flat_potrad[dayperyear] <- flat_potrad[dayperyear - 1]
  daylength[dayperyear] <- daylength[dayperyear - 1]
  tiny_rad_fract[dayperyear] <- tiny_rad_fract[dayperyear - 1]
  # return tiny_rad_fract, daylength, flat_potrad, tt_max0
  all<-NULL
  all$tt_max0<-tt_max0
  all$flat_potrad<-flat_potrad
  all$daylength<-daylength
  all$tiny_rad_fract<-tiny_rad_fract
  return(all)
}

solar_geom2 <- function(lat = lat) {
  # """
  # Flat earth assumption
  #
  # Parameters
  # ----------
  # lat:
  #     Latitude in decimal format
  #
  # Returns
  # -------
  # sg:
  #     (tiny_rad_fract, daylength, flat_potrad, tt_max0)
  # """
  
  # optical airmass by degrees
  OPTAM <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
             4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
             8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
  dayperyear <- ceiling(cnst$DAYS_PER_YEAR)
  tt_max0 <- rep(0, dayperyear)
  daylength <-  rep(0, dayperyear)
  flat_potrad <-  rep(0, dayperyear)
  
  # Translate lat to rad
  lat <- min(max(lat * cnst$RAD_PER_DEG, -pi / 2.), pi / 2.0)
  coslat <- cos(lat)
  sinlat <- sin(lat)
  
  # Sub-daily time step and angular step
  dt <- cnst$SW_RAD_DT
  dh <- dt / cnst$SEC_PER_RAD
  
  # Allocate the radiation arrays
  tiny_step_per_day <- cnst$SEC_PER_DAY / cnst$SW_RAD_DT
  tiny_rad_fract <- array(0, dim = c(dayperyear, tiny_step_per_day))
  for (i in 1:dayperyear) {
    # Declination and quantities of interest
    decl <- cnst$MIN_DECL * cos(((i-1) + cnst$DAYS_OFF) * cnst$RAD_PER_DAY)
    cosdecl <- cos(decl)
    sindecl <- sin(decl)
    
    # calculate daylength as a function of lat and decl
    cosegeom <- coslat * cosdecl
    sinegeom <- sinlat * sindecl
    coshss <- min(max(-sinegeom / cosegeom, -1), 1)
    hss <- acos(coshss)
    daylength[i] <- min(2.0 * hss * cnst$SEC_PER_RAD, cnst$SEC_PER_DAY)
    # Extraterrestrial radiation perpendicular to beam,
    # total over the timestep (J)
    dir_beam_topa <- (1368.0 + 45.5 * sin((2.0 * pi * (i-1) / cnst$DAYS_PER_YEAR) + 1.7)) * dt
    sum_trans <- 0
    sum_flat_potrad <- 0
    # Set up angular calculations
    for (h in seq(from = -hss, to = hss, by = dh)) {
      # Cosine of the hour angle and solar zenith angle
      cosh <- cos(h)
      cza <- cosegeom * cosh + sinegeom
      if (cza > 0) {
        # When sun is above flat horizon do flat-surface
        # calculations to determine daily total transmittance
        # and save potential radiation for calculation of
        # diffuse portion
        dir_flat_topa <- dir_beam_topa * cza
        am <- 1.0 / (cza + 0.0000001)
        if (am > 2.9) {
          ami <- min(max((acos(cza) / cnst$RAD_PER_DEG) - 69, 0), 20)
          # print(ami)
          am <- OPTAM[ami+1]
          # print(am)
        }
        sum_flat_potrad <- sum_flat_potrad + dir_flat_topa
      } else {
        # Sun not above horizon
        dir_flat_topa <- 0
      }
      tinystep <- min(max((12 * cnst$SEC_PER_HOUR + h * cnst$SEC_PER_RAD) / dt, 0),tiny_step_per_day - 1)
      tiny_rad_fract[i,tinystep] <- dir_flat_topa
    }
    if (daylength[i] && sum_flat_potrad > 0) {
      tiny_rad_fract[i,] <- tiny_rad_fract[i,] / sum_flat_potrad
    }
    if (daylength[i]) {
      # Transmittance and potential radiation
      # averaged over daylength
      flat_potrad[i] <- sum_flat_potrad / daylength[i]
    } else {
      # No daytime - no radiation
      tt_max0[i] <- 0.
      flat_potrad[i] <- 0.
    }
  }
  flat_potrad[dayperyear] <- flat_potrad[dayperyear - 1]
  daylength[dayperyear] <- daylength[dayperyear - 1]
  tiny_rad_fract[dayperyear] <- tiny_rad_fract[dayperyear - 1]
  # return tiny_rad_fract, daylength, flat_potrad, tt_max0
  all<-NULL
  all$flat_potrad<-flat_potrad
  all$daylength<-daylength
  all$tiny_rad_fract<-tiny_rad_fract
  return(all)
}

solar_geom_tt_max0 <- function(elev = elev, lat = lat, lr = lr) {
  # """
  # Flat earth assumption
  #
  # Parameters
  # ----------
  # elev:
  #     Elevation in meters
  # lat:
  #     Latitude in decimal format
  # lr:
  #     Lapse rate in K/m
  #
  # Returns
  # -------
  # sg:
  #     (tiny_rad_fract, daylength, flat_potrad, tt_max0)
  # """
  
  # optical airmass by degrees
  OPTAM <- c(2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07,
             4.37, 4.72, 5.12, 5.60, 6.18, 6.88, 7.77,
             8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00)
  dayperyear <- ceiling(cnst$DAYS_PER_YEAR)
  tt_max0 <- rep(0, dayperyear)
  daylength <-  rep(0, dayperyear)
  flat_potrad <-  rep(0, dayperyear)
  
  # Calculate pressure ratio as a function of elevation
  t1 <- 1.0 - (lr * elev) / cnst$T_STD
  t2 <- cnst$G_STD / (lr * (cnst$R / cnst$MA))
  trans <- cnst$TBASE^(t1^t2)
  
  # Translate lat to rad
  lat <- min(max(lat * cnst$RAD_PER_DEG, -pi / 2.), pi / 2.0)
  coslat <- cos(lat)
  sinlat <- sin(lat)
  
  # Sub-daily time step and angular step
  dt <- cnst$SW_RAD_DT
  dh <- dt / cnst$SEC_PER_RAD
  
  # Allocate the radiation arrays
  for (i in 1:dayperyear) {
    # Declination and quantities of interest
    decl <- cnst$MIN_DECL * cos(((i-1) + cnst$DAYS_OFF) * cnst$RAD_PER_DAY)
    cosdecl <- cos(decl)
    sindecl <- sin(decl)
    
    # calculate daylength as a function of lat and decl
    cosegeom <- coslat * cosdecl
    sinegeom <- sinlat * sindecl
    coshss <- min(max(-sinegeom / cosegeom, -1), 1)
    hss <- acos(coshss)
    daylength[i] <- min(2.0 * hss * cnst$SEC_PER_RAD, cnst$SEC_PER_DAY)
    # Extraterrestrial radiation perpendicular to beam,
    # total over the timestep (J)
    dir_beam_topa <- (1368.0 + 45.5 * sin((2.0 * pi * (i-1) / cnst$DAYS_PER_YEAR) + 1.7)) * dt
    sum_trans <- 0
    sum_flat_potrad <- 0
    # Set up angular calculations
    for (h in seq(from = -hss, to = hss, by = dh)) {
      # Cosine of the hour angle and solar zenith angle
      cosh <- cos(h)
      cza <- cosegeom * cosh + sinegeom
      if (cza > 0) {
        # When sun is above flat horizon do flat-surface
        # calculations to determine daily total transmittance
        # and save potential radiation for calculation of
        # diffuse portion
        dir_flat_topa <- dir_beam_topa * cza
        am <- 1.0 / (cza + 0.0000001)
        if (am > 2.9) {
          ami <- min(max((acos(cza) / cnst$RAD_PER_DEG) - 69, 0), 20)
          # print(ami)
          am <- OPTAM[ami+1]
          # print(am)
        }
        sum_trans <- sum_trans + ((trans^am) * dir_flat_topa)
        sum_flat_potrad <- sum_flat_potrad + dir_flat_topa
      } else {
        # Sun not above horizon
        dir_flat_topa <- 0
      }
    }
    if (daylength[i]) {
      # Transmittance and potential radiation
      # averaged over daylength
      tt_max0[i] <- sum_trans / sum_flat_potrad
    } else {
      # No daytime - no radiation
      tt_max0[i] <- 0.
    }
  }
  tt_max0[dayperyear] <- tt_max0[dayperyear - 1]
  return(tt_max0)
}

for (q in c(1:10)) {
  result <- solar_geom(900, 51.25, 0.0065)
}
csolar_geom <- cmpfun(solar_geom)

n <- 10
system.time(replicate(n,solar_geom(900, 51.25, 0.0065)))
system.time(replicate(n,csolar_geom(900, 51.25, 0.0065)))

csolar_geom_tt_max0 <- cmpfun(solar_geom_tt_max0)
system.time(replicate(n,solar_geom_tt_max0(900, 51.25, 0.0065)))
system.time(replicate(n,csolar_geom_tt_max0(900, 51.25, 0.0065)))

all1 <- solar_geom2(50.25)
all2 <- solar_geom2(50.25)
tt_max0 <- solar_geom_tt_max0(900, 51.25, 0.0065)
plot(result$tiny_rad_fract[70,])

# print(result$flat_potrad-all2$flat_potrad)
# print(result$daylength-all2$daylength)
# print(result$tiny_rad_fract[1,]-all2$tiny_rad_fract[1,])
# print(result$tiny_rad_fract[,100]-all2$tiny_rad_fract[,100])

print(result$tt_max0-tt_max0)


all2$daylength/3600

# arraysize_flat_potrad    <- array(0, dim = c(366, 360))
# arraysize_daylength      <- array(0, dim = c(366, 360))
# arraysize_tt_max0        <- array(0, dim = c(366, 360, 720)) # vanwege dat elke gridcell een andere elevation kan hebbben
# arraysize_tiny_rad_fract <- array(0, dim = c(366, 2880, 360))

for (q in c(1:10)) {
  result <- solar_geom(900, 51.25, 0.0065)
}
#################################################################################
library(WFRTools)
library(lattice)

outstep <- 6

## LOAD Precip
inElev <- ncLoad(file = "./WFDEI-elevation.nc",
                 varName = "elevation")

## output dimensions
stepspday <- 24/outstep # Define steps per day
nx <- length(inElev$xyCoords$x)
ny <- length(inElev$xyCoords$y)

## Define output arrays
tiny_step_per_day <- cnst$SEC_PER_DAY / cnst$SW_RAD_DT
outRadFractions <- array(data = NA, dim = c(366, tiny_step_per_day, ny))
outFlat_potrad  <- array(data = NA, dim = c(366, ny))
outDaylength    <- array(data = NA, dim = c(366, ny))
outTt_max0      <- array(data = NA, dim = c(366, ny, nx)) # vanwege dat elke gridcell een andere elevation kan hebbben

library(doParallel)
registerDoParallel(cores=3)

## Do solar_geom
# for (iy in (1:ny) ) {
#   for (ix in (1:nx) ) {
print(Sys.time())
for (iy in (1:1) ) {
  print(Sys.time())
  print(iy)
  # for (ix in (1:nx) ) {
  output<-foreach(ix = 1:10) %dopar% {
    if (!is.na(inElev$Data[iy,ix])) {
      lat <- inElev$xyCoords$y[iy]
      result <- solar_geom(inElev$Data[iy,ix], lat, 0.0065)
      # result <-ix
      print(result)
      # outRadFractions[, , iy] <- result$tiny_rad_fract
      # outFlat_potrad[, iy] <- result$flat_potrad
      # outDaylength[, iy] <- result$daylength
      # outTt_max0[, iy, ix] <- result$tt_max0
    }
  }
  for (ix in 1:nx) {
    if (!is.na(inElev$Data[iy,ix])) {
      ## ADD TO OUTPUT ARRAY
      outRadFractions[, , iy] <- output[[ix]]$tiny_rad_fract
      outFlat_potrad[, iy] <- output[[ix]]$flat_potrad
      outDaylength[, iy] <- output[[ix]]$daylength
      outTt_max0[, iy, ix] <- output[[ix]]$tt_max0
        # toNetCDFData[[iVar]][ix,iy,] <- output[[ix]][iStart:iEnd]
      # }
    }
  }

}
save(file = "../disaggregateR/outRadFractions.Rdata", outRadFractions)
save(file = "../disaggregateR/outFlat_potrad.Rdata", outFlat_potrad)
save(file = "../disaggregateR/outDaylength.Rdata", outDaylength)
save(file = "../disaggregateR/outTt_max0.Rdata", outTt_max0)




#
# ## Display for debugging
# levelplot(t(outRsds[1,,]), col.regions = rainbow(60))
# # ggplotje(inRsds)
#
# csolar_geom <- cmpfun(solar_geom)
# result <- solar_geom(900, 51.25, 0.0065)
