makeNetcdfOut <- function() {
  
  settings<-metGen$settings
  derived<-metGen$derived
  
  # CREATE NETCDF
  FillValue <- 1e20
  
  ## Define dimensions
  dimX <- ncdim_def("lon", "degrees_east", metGen$settings$x)
  dimY <- ncdim_def("lat", "degrees_north",metGen$settings$y)
  timeString <-format(strptime(settings$startDate, format = "%Y-%m-%d", tz = "GMT"),format="%Y-%m-%d %T")
  timeArray <-c(0:(metGen$derived$nrec_out-1)) * (24 / (24/metGen$derived$outDt))
  dimT <- ncdim_def("time", paste0("hours since ",timeString), timeArray, unlim = TRUE, calendar = "standard")
  
  dimsizes<-c(metGen$settings$nx,metGen$settings$ny,metGen$derived$nrec_out)
  ################
  for (var in names(settings$outVars)) {
    ## Create folder
    dir.create(file.path(getwd(), dirname(settings$outVars[[var]]$filename)), showWarnings = FALSE, recursive = T)
    
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

mgcheckInVars <- function() {
  for (var in names(metGen$settings$inVar)) {
    printf("Checking \"%s\"... ", var)
    ## Open the netcdf file
    ncFile <- nc_open( metGen$settings$inVar[[var]]$filename )
    attTmp<-ncatt_get( ncFile, metGen$settings$inVar[[var]]$ncname, "units" )
    if (attTmp$hasatt == TRUE) {
      metGen$settings$inVar[[var]]$units <- attTmp$value
    } else {
      metGen$settings$inVar[[var]]$units <- "missing"
    }
    
    ## Close the file
    nc_close(ncFile)
    
    unitIn<-metGen$settings$inVar[[var]]$units
    ## Correction of some commenly used units
    if (unitIn == "C") {
      if (verbose) printf("Unit is \"%s\" assuming that the unit is \"Celsius\". ", unitIn)
      unitIn <- "Celsius"
    }
    if (grepl("kg/m2", unitIn)) {
      unitInTmp<-gsub("kg/m2", "mm" ,unitIn)
      if (verbose) printf("Unit contains \"%s\" assuming that the unit is \"%s\". ", unitIn,unitInTmp)
      unitIn<-unitInTmp
    }
    if (grepl("kg m-2", unitIn)) {
      unitInTmp<-gsub("kg m-2", "mm" ,unitIn)
      if (verbose) printf("Unit contains \"%s\" assuming that the unit is \"%s\". ", unitIn,unitInTmp)
      unitIn<-unitInTmp
    }
    
    metGen$settings$inVar[[var]]$units<-unitIn
    convertUnit(data = NULL, metGen$settings$inVar[[var]]$units, metGen$metadata$invars[[var]]$units, verbose = T, doConversion = F)
    printf("\n")
  }
}

readAllForcing <- function(date) {
  ## Read data
  forcing_dataR <- NULL
  for (var in names(metGen$settings$inVar)) {

    forcing_dataR[[var]] <- ncLoad(filename = metGen$settings$inVar[[var]]$filename,
                                   var = metGen$settings$inVar[[var]]$ncname,
                                   lonlatbox = metGen$settings$lonlatbox,
                                   date = date)

    forcing_dataR[[var]] <- convertUnit(forcing_dataR[[var]],
                                        metGen$settings$inVar[[var]]$units,
                                        metGen$metadata$invars[[var]]$units)
  }
  
  return(forcing_dataR)
}


convertUnit <-function(data, unitIn, unitOut, verbose = F, doConversion = T) {
  if (!unitIn == unitOut) {
    if (ud.are.convertible(unitIn,unitOut)) {
      if (verbose) printf("\"%s\" will be converted to \"%s\". ", unitIn, unitOut)
      if (doConversion) data[]<-ud.convert(data[],unitIn,unitOut)
    } else {
      if (verbose) printf("\"%s\" cannot be converted to \"%s\". ", unitIn, unitOut)
    }
  } else {
    if (verbose) printf("\"%s\" does not need to be converted. ", unitIn)
  }
  if (doConversion) return(data)
}

##https://cran.r-project.org/web/packages/futureheatwaves/vignettes/starting_from_netcdf.html
ncLoad <- function(filename, var, lonlatbox, date = NULL) {
  
  lon_range <- c(lonlatbox[1], lonlatbox[2])
  lat_range <- c(lonlatbox[3], lonlatbox[4])
  
  ncid <- nc_open(filename = filename)
  lons <- ncvar_get(ncid, "lon")
  lats <- ncvar_get(ncid, "lat")
  
  # If requested, limit to certain ranges of latitude and/or longitude
  if(!is.null(lon_range)){
    lon_range <- sort(lon_range)
    lon_index <- which(lon_range[1] <= lons & lons <= lon_range[2]) 
    # lon <- lons[lon_index]
    # tas <- tas[lon_index, , ]
  }
  if(!is.null(lat_range)){
    lat_range <- sort(lat_range)
    lat_index <- which(lat_range[1] <= lats &  lats <= lat_range[2]) 
    # lat <- lats[lat_index]
    # tas <- tas[ , lat_index, ]
  }
  
  if (!is.null(date)) {
    times <- nc.get.time.series(ncid)
    time_index <- which(format(times, "%Y-%m-%d") == format(date, "%Y-%m-%d"))
    
    dataset <- nc.get.var.subset.by.axes(ncid, var,
                                         axis.indices = list(X = lon_index, 
                                                             Y = lat_index,
                                                             T = time_index)
                                         # axes.map = c(3,2,1)
    )
  } else {
    dataset <- nc.get.var.subset.by.axes(ncid, var,
                                         axis.indices = list(X = lon_index, 
                                                             Y = lat_index,
                                                             T = 1))
  }
  nc_close(ncid)
  
  ## Flip if needed
  if (lats[2] < lats[1])
    dataset[]<-dataset[,c(dim(dataset)[2]:1),]
  
  return(dataset)
}
