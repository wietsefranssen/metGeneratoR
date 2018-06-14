## Cleanup
rm(list=ls(all=TRUE))

library(metGeneratoR)

params <- NULL
params$time_step <- 60*3
params$do_rad_small <- F

cnst$SW_RAD_DT <- 30
# cnst$SW_RAD_DT <- 30*10

######### Do load rad_small
if(params$do_rad_small) {
  load(file = "./hpc/outDaylength_2880.Rdata")
  load(file = "./hpc/outRadFractions_2880.Rdata")
  
  rad_small <- rad_fract_chunk_small(outDaylength, outRadFractions, params)
  save(rad_small , file = "rad_small.Rdata")
  rm(outRadFractions, outDaylength)
} else {
  load(file = "rad_small.Rdata")
  load(file = "./hpc/outDaylength_2880.Rdata")
}

#############################################

## set the settings
settings <- mtclim_getSettings()
settings <- setOutstep(3, settings)
# settings$lonlatbox <- c(92.25, 110.25, 7.25, 36.25)
settings$lonlatbox <- c(-179.75, 179.75, -89.75, 89.75)
settings$system$nCores <- 4

settings <- setInputVars(settings,list(
  # tasmax         = list(ncFileName = "../metGeneratoR/hpc/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc",        ncName = "rsds", vicIndex = 16),
  # tas         = list(ncFileName = "../metGeneratoR/hpc/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc",        ncName = "rsds", vicIndex = 16),
  pr         = list(ncFileName = "../metGeneratoR/hpc/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc",        ncName = "rsds", vicIndex = 16)
))
settings$elevation <- list(ncFileName = "./WFDEI-elevation.nc", ncName = "elevation")

## Output variables
## Comment out the ones you dont want to include
settings <- setOutputVar(settings, "shortwave")
settings <- setOutputVar(settings, "pr")

# ## Comment out the ones you dont want to include
# settings$outputVars$shortwave$do <- list(
#   # pr         = list(VICName = "OUT_PREC",       units = "mm",        longName = "incoming precipitation"),
#   # tas        = list(VICName = "OUT_AIR_TEMP",   units = "C",         longName = "air temperature"),
#   shortwave  = list(VICName = "OUT_SHORTWAVE",  units = "W m-2",     longName = "incoming shortwave"),
#   # longwave   = list(VICName = "OUT_LONGWAVE",   units = "W m-2",     longName = "incoming longwave"),
#   # pressure   = list(VICName = "OUT_PRESSURE",   units = "kPa",       longName = "near surface atmospheric pressure"),
#   # qair       = list(VICName = "OUT_QAIR",       units = "kg kg-1",   longName = "specific humidity"),
#   # vp         = list(VICName = "OUT_VP",         units = "kPa",       longName = "near surface vapor pressure"),
#   # rel_humid  = list(VICName = "OUT_REL_HUMID",  units = "fraction",  longName = "relative humidity"),
#   # density    = list(VICName = "OUT_DENSITY",    units = "kg m-3",    longName = "near-surface atmospheric density"),
#   wind       = list(VICName = "OUT_WIND",       units = "m s-1",     longName = "near surface wind speed")
# )

## Start profiler
start.time.total <- Sys.time()

## Register nr of cores
cat(paste("nCores: ", settings$system$nCores),"\n")
registerDoParallel(cores=settings$system$nCores)

## Set outvars in settings
settings$mtclim$nOut <- length(settings$outputVars)
for (i in 1:length(settings$outputVars)) {
  settings$mtclim$outNames[i]<-settings$outputVars[[i]]$VICName
}

## LOAD MASK/ELEVATION
elevation <- ncLoad(file = settings$elevation$ncFileName,
                    varName = settings$elevation$ncName,
                    lonlatbox = settings$lonlatbox)
mask<-elevation

## Subselect from global data
startX <- ((settings$lonlatbox[1] - -179.75) * 2) + 1
endX <- ((settings$lonlatbox[2] - -179.75) * 2) + 1
startY <- ((settings$lonlatbox[3] - -89.75) * 2) + 1
endY <- ((settings$lonlatbox[4] - -89.75) * 2) + 1
rad_small <- rad_small[startY:endY,,]
outDaylength <- outDaylength[startY:endY,]

## makeOutputNetCDF
for (iYear in 1:length(settings$ncOut)) {
  makeNetcdfOut(settings, elevation, settings$ncOut[[iYear]], "shortwave")
  makeNetcdfOut(settings, elevation, settings$ncOut[[iYear]], "pr")
}

## Change settings for current part
elevation <- ncLoad(file = settings$elevation$ncFileName,
                    varName = settings$elevation$ncName,
                    lonlatbox = settings$lonlatbox)

nx <- length(elevation$xyCoords$x)
ny <- length(elevation$xyCoords$y)

## Print part info
cat(sprintf("\nRunning\n"))

## DEFINE OUTPUT ARRAY
el <- array(NA, dim = c(nx, ny, settings$outstep_per_day))
toNetCDFData <- list(el)[rep(1,length(settings$outputVars))]
outShortwave <- el
rm(el)

rad_small <- rad_small * settings$outstep_per_day

profile<-NULL
for (iday in 1:settings$intern$nrec_in) {
  ## LOAD WHOLE DOMAIN FROM NETCDF
  profile$start.time.read <- Sys.time()
  forcing_dataRTotal <- readAllForcing(settings, elevation, iday)
  profile$end.time.read <- Sys.time()
  
  ## Init progressbar
  pb <- txtProgressBar(min = 0, max = ny, initial = 0, char = ">",
                       width = 80, title, label, style = 1, file = "")
  
  profile$start.time.run <- Sys.time()
  
  ## Cut rad_small per day
  # if (settings$outputVars$shortwave)
  rad_small_xy <- array(NA, dim = c(nx, ny, settings$outstep_per_day))  
  for (ix in 1:nx) {
    rad_small_xy[ix,,] <- rad_small[ , iday, ]
  }
  
  ## Do radiation
  if (settings$outputVars$shortwave$enable) {
    for (iVar in 1:length(settings$outputVars)) {
      for (it in 1:settings$outstep_per_day) {
        outShortwave[,,it] <- forcing_dataRTotal[[1]][ , , 1 ] * rad_small_xy[,,it] 
      }
    }
  }
  profile$end.time.run <- Sys.time()
  
  # Close ProgressBar
  close(pb)
  
  ## ADD OUTPUT TO NETCDF
  ## Define numer of years
  profile$start.time.write <- Sys.time()
  for (iYear in 1:length(settings$ncOut)) {
    if (settings$outputVars$shortwave$enable) {
      timeIndex <- settings$outstep_per_day*(iday-1)+1
      ncid <- nc_open(settings$ncOut[[iYear]]$fileName, write = TRUE)
      ncvar_put(ncid,
                "shortwave",
                outShortwave[,,],
                start = c(1, 1, timeIndex),
                count = c(nx, ny, settings$outstep_per_day)
      )
      nc_close(ncid)
    }
    if (settings$outputVars$pr$enable) {
      timeIndex <- settings$outstep_per_day*(iday-1)+1
      ncid <- nc_open(settings$ncOut[[iYear]]$fileName, write = TRUE)
      ncvar_put(ncid,
                "pr",
                outPr[,,],
                start = c(1, 1, timeIndex),
                count = c(nx, ny, settings$outstep_per_day)
      )
      nc_close(ncid)
    }
  }
  
  # timeString <-format(strptime(ncOut$startdate, format = "%Y-%m-%d", tz = "GMT"),format="%Y-%m-%d %T")
  # timeArray <-c(0:(ncOut$nrec_out-1)) * (24 / (24/settings$outstep))
  # dimT <- ncdim_def("time", paste0("hours since ",timeString), 0, unlim = TRUE, calendar = "standard")
  
  # for (iYear in 1:length(settings$ncOut)) {
  #   ncid <- nc_open(settings$ncOut[[iYear]]$fileName, write = TRUE)
  #   for (iVar in 1:length(settings$outputVars))
  #   {
  #     timeIndex <- settings$outstep_per_day*(iday-1)+1
  #     ncvar_put(ncid,
  #               names(settings$outputVars)[iVar],
  #               toNetCDFData[[iVar]][,,],
  #               start = c(1,
  #                         1,
  #                         timeIndex),
  #               count = c(nx,
  #                         ny,
  #                         settings$outstep_per_day)
  #     )
  #   }
  #   nc_close(ncid)
  # }
  profile$end.time.write <- Sys.time()
  
  ## Print info about part
  cat(sprintf("  Times (read/run/write/total): %.1f/%.1f/%.1f/%.1f seconds\n",
              as.numeric(profile$end.time.read  - profile$start.time.read, units = "secs"),
              as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs"),
              as.numeric(profile$end.time.write - profile$start.time.write, units = "secs"),
              as.numeric(profile$end.time.write - profile$start.time.read, units = "secs")))
  cat(sprintf("  Sizes (read/write): %s/%s\n",
              format(object.size(forcing_dataRTotal), units = "auto"),
              format(object.size(toNetCDFData), units = "auto")))
}

cat(sprintf("\nFinished in %.1f minutes\n",as.numeric(Sys.time() - start.time.total, units = "mins")))
