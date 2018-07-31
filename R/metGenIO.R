makeNetcdfOut <- function(mask) {
  
  settings<-metGen$settings
  derived<-metGen$derived
  
  # CREATE NETCDF
  FillValue <- 1e20
  
  ## Define dimensions
  dimX <- ncdim_def("lon", "degrees_east", mask$xyCoords$x)
  dimY <- ncdim_def("lat", "degrees_north",mask$xyCoords$y)
  timeString <-format(strptime(settings$startDate, format = "%Y-%m-%d", tz = "GMT"),format="%Y-%m-%d %T")
  timeArray <-c(0:(metGen$derived$nrec_out-1)) * (24 / (24/metGen$derived$outDt))
  dimT <- ncdim_def("time", paste0("hours since ",timeString), timeArray, unlim = TRUE, calendar = "standard")
  
  dimsizes<-c(length(mask$xyCoords$x),length(mask$xyCoords$y),metGen$derived$nrec_out)
  ################
  for (var in names(settings$outVars)) {
    ## Create folder
    dir.create(file.path(getwd(), basename(dirname(settings$outVars[[var]]$filename))), showWarnings = FALSE)
    
    dataVar <- ncvar_def(name=var, units='', compression = 7, dim=list(dimX,dimY,dimT), missval=FillValue, prec="float")
  
    ## SAVE AS NC-DATA
    cat(sprintf("Create output file: %s\n", settings$outVars[[var]]$filename))
    ncid <- nc_create(settings$outVars[[var]]$filename, dataVar, force_v4=TRUE)
    
    ncatt_put( ncid, "lon", "standard_name", "longitude")
    ncatt_put( ncid, "lon", "long_name",     "Longitude")
    ncatt_put( ncid, "lon", "axis",          "X")
    ncatt_put( ncid, "lat", "standard_name", "latitude")
    ncatt_put( ncid, "lat", "long_name",     "Latitude")
    ncatt_put( ncid, "lat", "axis",          "Y")
    ncatt_put( ncid, "time", "standard_name", "time")
    ncatt_put( ncid, "time", "calendar",     "standard")
    ncatt_put( ncid, names(settings$outVars[var]), "standard_name", names(settings$outVars[var]))
    ncatt_put( ncid, names(settings$outVars[var]), "long_name", settings$outVars[[var]]$longName)
    ncatt_put( ncid, names(settings$outVars[var]), "units", settings$outVars[[var]]$units)
    
    # ## Global Attributes
    ncatt_put( ncid, 0, "NetcdfCreatationDate", as.character(Sys.Date()))
    
    ## Close Netcdf file
    nc_close(ncid)
  }
}

readAllForcing <- function(mask, timestep) {
  # mask<-elevation
  nx <- length(mask$xyCoords$x)
  ny <- length(mask$xyCoords$y)
  
  # forcing_dataR <- list()
  # for (i in 1:length(metGen$settings$inVar)) {
  #   forcing_dataR[[i]]<- array(0, dim=c(nx, ny, 1))
  # }
  
  ## Read data
  forcing_dataR <- NULL
  for (var in names(metGen$settings$inVar)) {
    # for (var in 1:length(metGen$settings$inVar)) {
    # forcing_dataR[[var]][] <- ncLoad(file = metGen$settings$inVar[[var]]$filename,
                                     forcing_dataR[[var]] <- ncLoad(file = metGen$settings$inVar[[var]]$filename,
                                                                      varName = metGen$settings$inVar[[var]]$ncname,
                                      # lonlatbox = c(settings$lonlatbox[1],
                                      #               settings$lonlatbox[2],
                                      #               settings$lonlatbox[3],
                                      #               settings$lonlatbox[4]),
                                      lonlatbox = metGen$settings$lonlatbox,
                                      timesteps = timestep)$Data
  }
  return(forcing_dataR)
}
