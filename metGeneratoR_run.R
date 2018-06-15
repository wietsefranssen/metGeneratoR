## Cleanup
rm(list=ls(all=TRUE))

library(metGeneratoR)

# source("./R/metGenSettings.R")

mgsetLonlatbox(c(92.25, 110.25, 7.25, 36.25))
# mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))
mgsetPeriod(startdate = "1950-01-01", enddate = "1950-01-02")
mgsetNHourPerStep(3) # Set N hours per timestep
mgsetNCores(4) # Set N hours per timestep

## Define input variables
# mgsetInVars(list(
#   pr         = list(ncname = "pr",     filename = "./inst/extdata/pr_Mekong.nc4"),
#   tasmax     = list(ncname = "tasmax", filename = "./inst/extdata/tasmax_Mekong.nc4"),
#   wind       = list(ncname = "wind", filename = "./inst/extdata/wind_Mekong.nc4")
# ))

mgsetInVars(list(
  pr         = list(ncname = "pr",     filename = "../example_data4mtclim/Global/pr_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  tasmin     = list(ncname = "tasmin", filename = "../example_data4mtclim/Global/tasmin_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  tasmax     = list(ncname = "tasmax", filename = "../example_data4mtclim/Global/tasmax_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  shortwave  = list(ncname = "rsds",   filename = "../example_data4mtclim/Global/rsds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  longwave   = list(ncname = "rlds",   filename = "../example_data4mtclim/Global/rlds_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc"),
  wind       = list(ncname = "sfcWind",   filename = "../example_data4mtclim/Global/wind_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_1998.nc")
))

## Define elevation file
mgsetElevation(ncname = "elevation", filename = metGen$internal$ncFileNameElevation)

## Define output variables
mgsetOutVars(c("pr", "tas"))
# mgsetOutVars(c( "shortwave"))



##############################################
## Start profiler
start.time.total <- Sys.time()

## Register nr of cores
cat(paste("nCores: ", metGen$settings$nCores),"\n")
registerDoParallel(cores=metGen$settings$nCores)

## LOAD MASK/ELEVATION
elevation <- ncLoad(file = metGen$internal$ncFileNameElevation,
                    varName = metGen$settings$elevation$ncname,
                    lonlatbox = metGen$settings$lonlatbox)
mask<-elevation
nx <- length(mask$xyCoords$x)
ny <- length(mask$xyCoords$y)

## Do load rad_small & Subselect from global data
do_rad_small <- F
if(do_rad_small) {
  load(file = "./hpc/outDaylength_2880.Rdata")
  load(file = "./hpc/outRadFractions_2880.Rdata")
  
  rad_small <- rad_fract_chunk_small(outDaylength, outRadFractions, params) ## TODO!!
  save(rad_small , file = "rad_small.Rdata")
  rm(outRadFractions, outDaylength)
} else {
  load(file = "rad_small.Rdata")
  load(file = "./hpc/outDaylength_2880.Rdata")
}

startX <- ((metGen$settings$lonlatbox[1] - -179.75) * 2) + 1
endX <- ((metGen$settings$lonlatbox[2] - -179.75) * 2) + 1
startY <- ((metGen$settings$lonlatbox[3] - -89.75) * 2) + 1
endY <- ((metGen$settings$lonlatbox[4] - -89.75) * 2) + 1
rad_small <- rad_small[startY:endY,,]
rad_small <- rad_small * metGen$settings$nOutStepDay
outDaylength <- outDaylength[startY:endY,]

## makeOutputNetCDF
makeNetcdfOut(mask)

## DEFINE OUTPUT ARRAY
outData <- NULL
for (var in names(metGen$settings$outVars)) {
  outData[[var]] <- array(NA, dim = c(nx, ny, metGen$settings$nOutStepDay))
}

profile<-NULL
for (iday in 1:metGen$derived$nrec_in) {
  cat(paste("Running timestep: ", iday))
  
  ## LOAD WHOLE DOMAIN FROM NETCDF
  profile$start.time.read <- Sys.time()
  inData <- readAllForcing(elevation, iday)
  profile$end.time.read <- Sys.time()
  
  ## Init progressbar
  pb <- txtProgressBar(min = 0, max = ny, initial = 0, char = ">",
                       width = 80, title, label, style = 1, file = "")
  
  profile$start.time.run <- Sys.time()
  
  ## Cut rad_small per day
  rad_small_xy <- array(NA, dim = c(nx, ny, metGen$settings$nOutStepDay))
  for (ix in 1:nx) {
    rad_small_xy[ix,,] <- rad_small[ , iday, ]
  }
  
  ## Do precip
  if (!is.null(metGen$settings$outVars$pr)) {
    for (it in 1:metGen$settings$nOutStepDay) {
      outData[["pr"]][,,it] <- inData[["pr"]][ , , 1 ]
    }
  }
  
  ## Do radiation
  if (!is.null(metGen$settings$outVars$shortwave)) {
    for (it in 1:metGen$settings$nOutStepDay) {
      outData[["shortwave"]][,,it] <- inData[[1]][ , , 1 ] * rad_small_xy[,,it]
    }
  }
  profile$end.time.run <- Sys.time()
  
  # Close ProgressBar
  close(pb)
  
  ## ADD OUTPUT TO NETCDF
  profile$start.time.write <- Sys.time()
  for (var in names(metGen$settings$outVars)) {
    timeIndex <- metGen$settings$nOutStepDay*(iday-1)+1
    metGen$settings$outVars[[var]]$ncid <- nc_open(metGen$settings$outVars[[var]]$filename, write = TRUE)
    ncvar_put(metGen$settings$outVars[[var]]$ncid,
              var,
              outData[[var]][,,],
              start = c(1, 1, timeIndex),
              count = c(nx, ny, metGen$settings$nOutStepDay)
    )
    nc_close(metGen$settings$outVars[[var]]$ncid)
  }
  profile$end.time.write <- Sys.time()
  
  ## Print info about part
  cat(sprintf("  Times (read/run/write/total): %.1f/%.1f/%.1f/%.1f seconds\n",
              as.numeric(profile$end.time.read  - profile$start.time.read, units = "secs"),
              as.numeric(profile$end.time.run   - profile$start.time.run, units = "secs"),
              as.numeric(profile$end.time.write - profile$start.time.write, units = "secs"),
              as.numeric(profile$end.time.write - profile$start.time.read, units = "secs")))
  cat(sprintf("  Sizes (read/write): %s/%s\n",
              format(object.size(inData), units = "auto"),
              format(object.size(outData), units = "auto")
              # format(object.size(toNetCDFData), units = "auto")
  ))
}

cat(sprintf("\nFinished in %.1f minutes\n",as.numeric(Sys.time() - start.time.total, units = "mins")))
