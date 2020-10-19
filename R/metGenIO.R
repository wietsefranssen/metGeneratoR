#' @export
mggetInVarInfo <- function(filename, varname) {
  getdimname <- function(id) {
    result <- NULL
    i <- 1 
    for (iid in c(id)) {
      lijst <- c(dim$t$id,dim$x$id,dim$y$id)
      resultid<-which(lijst == iid)
      result[i] <- dim[[resultid]]$name
      i <- i + 1
    }
    return(result)
  }
  
  dimnames <- list(t="time", x="lon", y="lat")
  dimvarnames <- list(t="time", x="lon", y="lat")
  datavarname <- varname
  
  ncid <- open.nc(filename)
  
  result <- NULL
  
  ## Get dims
  dim <- NULL
  for (idim in c(1:length(dimnames))) {
    info <- dim.inq.nc(ncid, dimnames[[idim]])
    n <- info$length
    id <- info$id
    name <- info$name
    dim[[idim]] <- list(n=n, id=id, name=name)
  }
  
  ## change/add names to list members
  names(dim) <-  names(dimnames)
  
  ## add to result
  result$dims <- dim
  
  ## get vars
  varnames <- c(dimvarnames, varname)
  var <- NULL
  for (i in c(1:length(varnames))) {
    info <- var.inq.nc(ncid, varnames[[i]])
    name <- info$name
    id <- info$id
    ndim <- info$ndims
    dimids <- info$dimids
    dimnames <- getdimname(dimids)
    natts <- info$natts
    if (i != length(varnames)) {
      vals <- var.get.nc(ncid, varnames[[i]])
      var[[i]] <- list(id=id, name=name, ndim=ndim, dimids=dimids, dimnames=dimnames, natts=natts, vals=vals)
    } else {
      vals <- var.get.nc(ncid, varnames[[i]])
      var[[i]] <- list(id=id, name=name, ndim=ndim, dimids=dimids, dimnames=dimnames, natts=natts, vals=vals)
    }
  }
  
  ## change/add names to list members
  names(var) <-  c(names(dimvarnames), "data")
  
  ## add to result
  result$vars <- var
  
  close.nc(ncid)
  
  return(result)
}

#' @export
mggetInDims <- function() {
  
  for (var in names(metGen$settings$inVar)) {
    if (var != "radfrac") {
      printf("Getting input dimensions \"%s\". ", var)
      metGen$input[[var]] <-mggetInVarInfo(metGen$settings$inVars[[var]]$filename, metGen$settings$inVars[[var]]$ncname)
      # ncid <- nc_open(filename = metGen$settings$inVars[[var]]$filename)
      # 
      # metGen$input[[var]]$lons <- ncvar_get(ncid, "lon")
      # metGen$input[[var]]$lats <- ncvar_get(ncid, "lat")
      # metGen$input[[var]]$ndim <- length(dim(metGen$input[[var]]$lons))
      # nc_close(ncid)
      # 
      printf("\n")
    }
  }
}

#' @export
makeNetcdfOut <- function() {
  
  settings<-metGen$settings
  derived<-metGen$derived
  
  # CREATE NETCDF
  FillValue <- 1e20
  
  ## Define dimensions
  # dimX <- ncdim_def("lon", "degrees_east", metGen$settings$x)
  # dimY <- ncdim_def("lat", "degrees_north",metGen$settings$y)
  ################
  for (var in names(settings$outVars)) {
    if (metGen$output$ndim == 2) {
      dimX <- ncdim_def("lon", units = '', vals = c(1:length(metGen$output$lons[1,])), unlim = F, create_dimvar = F)
      dimY <- ncdim_def("lat", units = '', vals = c(1:length(metGen$output$lats[,1])), unlim = F, create_dimvar = F)
    } else {
      dimX <- ncdim_def("lon", units = '', vals = c(1:length(metGen$output$lons)), unlim = F, create_dimvar = F)
      dimY <- ncdim_def("lat", units = '', vals = c(1:length(metGen$output$lats)), unlim = F, create_dimvar = F)
    }
    timeString <-format(strptime(settings$startDate, format = "%Y-%m-%d", tz = "GMT"),format="%Y-%m-%d %T")
    timeArray <-c(0:(metGen$derived$nrec_out-1)) * (24 / (24/metGen$derived$outDt))
    dimT <- ncdim_def("time", paste0("hours since ",timeString), timeArray, unlim = TRUE, calendar = "standard")
    
    dimsizes<-c(metGen$settings$nx,metGen$settings$ny,metGen$derived$nrec_out)
    ## Create folder
    dir.create(file.path(getwd(), dirname(settings$outVars[[var]]$filename)), showWarnings = FALSE, recursive = T)
    
    vars <- NULL
    if (metGen$output$ndim == 2) {
      vars$lonVar <- ncvar_def(name='lon', units='', compression = 7, dim=list(dimY,dimX), missval=FillValue, prec="float")
      vars$latVar <- ncvar_def(name='lat', units='', compression = 7, dim=list(dimY,dimX), missval=FillValue, prec="float")
    } else {
      vars$lonVar <- ncvar_def(name='lon', units='', compression = 7, dim=dimX, missval=FillValue, prec="float")
      vars$latVar <- ncvar_def(name='lat', units='', compression = 7, dim=dimY, missval=FillValue, prec="float")
    }      
    vars$dataVar <- ncvar_def(name=var, units='', compression = 7, dim=list(dimY,dimX,dimT), missval=FillValue, prec="float")
    
    ## SAVE AS NC-DATA
    cat(sprintf("Create output file: %s\n", settings$outVars[[var]]$filename))
    ncid <- nc_create(settings$outVars[[var]]$filename, vars = vars, force_v4=TRUE)
    
    ncvar_put(ncid, 'lon', metGen$output$lons)
    ncvar_put(ncid, 'lat', metGen$output$lats)
    
    ncatt_put( ncid, "lon", "standard_name", "longitude")
    ncatt_put( ncid, "lon", "long_name",     "longitude")
    ncatt_put( ncid, "lon", "units",          "degrees_east")
    ncatt_put( ncid, "lon", "_CoordinateAxisType",          "Lon")
    
    ncatt_put( ncid, "lat", "standard_name", "latitude")
    ncatt_put( ncid, "lat", "long_name",     "Latitude")
    ncatt_put( ncid, "lat", "units",          "degrees_north")
    ncatt_put( ncid, "lat", "_CoordinateAxisType",          "Lat")
    
    ncatt_put( ncid, "time", "standard_name", "time")
    ncatt_put( ncid, "time", "calendar",     "standard")
    ncatt_put( ncid, names(settings$outVars[var]), "standard_name", names(settings$outVars[var]))
    ncatt_put( ncid, names(settings$outVars[var]), "long_name", settings$outVars[[var]]$longName)
    ncatt_put( ncid, names(settings$outVars[var]), "units", settings$outVars[[var]]$units)
    
    ## Global Attributes
    ncatt_put( ncid, 0, "ORIGINAL" , "*** ORIGINAL ATTRIBUTES ***")
    ncidIn<-nc_open(metGen$settings$inVars[[1]]$filename)
    nc.copy.atts(ncidIn, 0, ncid, 0)
    nc_close(ncidIn)
    ncatt_put( ncid, 0, "MetGenerator_package" , "This file is generated by the r-package: metGenerator")
    ncatt_put( ncid, 0, "MetGenerator_contact" , "Wietse Franssen (wietse.franssen@wur.nl)")
    ncatt_put( ncid, 0, "MetGenerator_created" , paste0(as.character(Sys.Date())))
    
    ## Close Netcdf file
    nc_close(ncid)
  }
}

#' @export
mgcheckInVars <- function() {
  for (var in names(metGen$settings$inVars)) {
    printf("Checking input variable \"%s\". ", var)
    ## Open the netcdf file
    ncFile <- nc_open( metGen$settings$inVars[[var]]$filename )
    attTmp <- ncatt_get( ncFile, metGen$settings$inVars[[var]]$ncname, "units" )
    if (attTmp$hasatt == TRUE) {
      metGen$metadata$inVars[[var]]$input_units <- attTmp$value
    } else {
      metGen$metadata$inVars[[var]]$input_units <- "missing"
    }
    
    ## Close the file
    nc_close(ncFile)
    
    unitIn<-metGen$metadata$inVars[[var]]$input_units
    
    ## Correction of some commenly used units
    if (unitIn == "C") {
      if (verbose) printf("Unit is \"%s\" assuming that the unit is \"Celsius\". ", unitIn)
      unitIn <- "Celsius"
    }
    if (grepl("kg/m2s", unitIn)) {
      unitInTmp <- gsub("kg/m2s", "mm s-1" ,unitIn)
      if (verbose) printf("Unit contains \"%s\" assuming that the unit is \"%s\". ", unitIn,unitInTmp)
      unitIn<-unitInTmp
    }
    if (grepl("kg/m2", unitIn)) {
      unitInTmp <- gsub("kg/m2", "mm" ,unitIn)
      if (verbose) printf("Unit contains \"%s\" assuming that the unit is \"%s\". ", unitIn,unitInTmp)
      unitIn<-unitInTmp
    }
    if (grepl("kg m-2", unitIn)) {
      unitInTmp <- gsub("kg m-2", "mm" ,unitIn)
      if (verbose) printf("Unit contains \"%s\" assuming that the unit is \"%s\". ", unitIn,unitInTmp)
      unitIn<-unitInTmp
    }
    
    metGen$metadata$inVars[[var]]$input_units <- unitIn
    
    convertUnit(data = NULL,
                metGen$metadata$inVars[[var]]$input_units,
                metGen$metadata$inVars[[var]]$internal_units,
                metGen$settings$inDt,
                verbose = T, doConversion = F)
    printf("\n")
  }
}

#' @export
mgcheckOutVars <- function() {
  for (var in names(metGen$settings$outVar)) {
    printf("Checking output variable \"%s\". ", var)
    
    convertUnit(outData[[var]], 
                metGen$metadata$outVars[[var]]$internal_units, 
                metGen$metadata$outVars[[var]]$output_units,
                metGen$settings$outDt,
                verbose = T, doConversion = F)
    
    printf("\n")
  }
}


#' @export
mgsetOutDims <- function() {
  printf("Setting output dimensions. ")
  
  if (is.null(metGen$settings$xybox[1])) {
    xybox <- NULL
    xybox[2] <- metGen$input[[1]]$dims$x$n
    xybox[4] <- metGen$input[[1]]$dims$y$n
    mgsetXYbox(c(1,xybox[2],1,xybox[4]))
  } 
  xx <- c(metGen$settings$xybox[1]:metGen$settings$xybox[2])
  yy <- c(metGen$settings$xybox[3]:metGen$settings$xybox[4])
  
  if (metGen$input[[1]]$vars$x$ndim == 1) {
    metGen$output$lons <- metGen$input[[1]]$vars$x$vals[xx]
    metGen$output$lats <- metGen$input[[1]]$vars$y$vals[yy] 
    metGen$output$ndim <- metGen$input[[1]]$vars$x$ndim
  } else if (metGen$input[[1]]$vars$x$ndim == 2) {
    if (metGen$input[[1]]$dims$x$id == metGen$input[[1]]$vars$x$dimids[1]) {
      metGen$output$lons <- metGen$input[[1]]$vars$x$vals[xx,yy]
      metGen$output$lats <- metGen$input[[1]]$vars$y$vals[xx,yy]
    } else if (metGen$input[[1]]$dims$x$id == metGen$input[[1]]$vars$x$dimids[2]) {
      metGen$output$lons <- metGen$input[[1]]$vars$x$vals[yy,xx]
      metGen$output$lats <- metGen$input[[1]]$vars$y$vals[yy,xx]
    } else {
      stop("error!")
    }
    metGen$output$ndim <- metGen$input[[1]]$vars$x$ndim
  } else {
    stop("too many dimensions!")
  }
  
  printf("\n")
  
  
}

#' @export
readAllForcing <- function(date) {
  ## Read data
  forcing_dataR <- NULL
  for (var in names(metGen$settings$inVars)) {
    forcing_dataR[[var]] <- ncLoad(filename = metGen$settings$inVars[[var]]$filename,
                                   var = metGen$settings$inVars[[var]]$ncname,
                                   xybox = metGen$settings$xybox,
                                   date = date)
    
    forcing_dataR[[var]] <- convertUnit(forcing_dataR[[var]],
                                        metGen$metadata$inVars[[var]]$input_units,
                                        metGen$metadata$inVars[[var]]$internal_units,
                                        metGen$settings$inDt)
  }
  
  return(forcing_dataR)
}

#' @export
convertUnit <- function(data, unitIn, unitOut, dt = 24, verbose = F, doConversion = T) {
  if (!unitIn == unitOut) {
    if (ud.are.convertible(unitIn, unitOut)) {
      if (verbose) printf("\"%s\" will be converted to \"%s\". ", unitIn, unitOut)
      if (doConversion) data[] <- ud.convert(data[],unitIn,unitOut)
    } else {
      if (unitIn == "mm s-1" && unitOut == "mm") {
        if (verbose) printf("\"%s\" will be converted to \"%s\". ", unitIn, unitOut)
        if (doConversion) {
          # print(data[5:10,312:315,1])
          # data[] <- ud.convert(data[],unitIn,unitOut)
          data[] <- data[] * (3600 * dt)
          # print(data[5:10,312:315,1])
        }
      } else if (unitIn == "mm" && unitOut == "mm s-1") {
        if (verbose) printf("\"%s\" will be converted to \"%s\". ", unitIn, unitOut)
        if (doConversion) {
          # print(data[5:10,312:315,1])
          # data[] <- ud.convert(data[],unitIn,unitOut)
          data[] <- data[] / (3600 * dt)
          # print(data[5:10,312:315,1])
        }
      } else {
        if (verbose) printf("\"%s\" cannot be converted to \"%s\". ", unitIn, unitOut)
      }
    }
  } else {
    if (verbose) printf("\"%s\" does not need to be converted. ", unitIn)
  }
  if (doConversion) return(data)
}

##https://cran.r-project.org/web/packages/futureheatwaves/vignettes/starting_from_netcdf.html
#' @export
ncLoad <- function(filename, var, xybox, date = NULL) {
  
  # lon_range <- c(xybox[1], xybox[2])
  # lat_range <- c(xybox[3], xybox[4])
  
  ncid <- nc_open(filename = filename)
  # lons <- ncvar_get(ncid, "lon")
  # lats <- ncvar_get(ncid, "lat")
  
  # # If requested, limit to certain ranges of latitude and/or longitude
  # if(!is.null(lon_range)){
  #   lon_range <- sort(lon_range)
  #   lon_index <- which(lon_range[1] <= lons & lons <= lon_range[2]) 
  #   # lon <- lons[lon_index]
  #   # tas <- tas[lon_index, , ]
  # }
  # if(!is.null(lat_range)){
  #   lat_range <- sort(lat_range)
  #   lat_index <- which(lat_range[1] <= lats &  lats <= lat_range[2]) 
  #   # lat <- lats[lat_index]
  #   # tas <- tas[ , lat_index, ]
  # }
  
  if (!is.null(date)) {
    ## TODO Move this part to a general part like check input data
    ## mgsetInDt(inDt = 3) this function can then also be filled automatically
    if (strsplit(ncid$dim$time$units,split=" ")[[1]][1] == "seconds") {
      datetmp <- strsplit(ncid$dim$time$units,split=" ")[[1]][3]
      timetmp <- strsplit(ncid$dim$time$units,split=" ")[[1]][4]
      times <- as.PCICt(paste(datetmp,timetmp), cal = "gregorian")+ncid$dim$time$vals
    } else {
      times <- nc.get.time.series(ncid)
    }
    
    time_index <- which(format(times, "%Y-%m-%d") == format(date, "%Y-%m-%d"))
    
    dataset <-  ncvar_get(nc = ncid, 
                          varid = var, 
                          start = c(xybox[3],xybox[1],time_index[1]),
                          count = c(length(c(xybox[3]:xybox[4])),
                                    length(c(xybox[1]:xybox[2])),length(time_index)),
                          collapse_degen = F)
    print(dim(dataset))
    # dataset <- nc.get.var.subset.by.axes(ncid, var,
    #                                      axis.indices = list(X = c(xybox[1]:xybox[2]), 
    #                                                          Y = c(xybox[3]:xybox[4]),
    #                                                          T = time_index,
    #                                      axes.map = c(0,1,2))
    # )
  } else {
    dataset <- nc.get.var.subset.by.axes(ncid, var,
                                         axis.indices = list(X = c(xybox[1]:xybox[2]), 
                                                             Y = c(xybox[3]:xybox[4]),
                                                             T = 1))
  }
  nc_close(ncid)
  
  # ## Flip if needed
  # if (lats[2] < lats[1])
  #   dataset[]<-dataset[,c(dim(dataset)[2]:1),]
  
  return(dataset)
}
