melo_create_nc <-function(nc_info) {
  ## Get nc_info
  nx <- nc_info$nx
  ny <- nc_info$ny
  lon_2d <- nc_info$lon_2d
  lat_2d <- nc_info$lat_2d
  inx <- nc_info$inx
  iny <- nc_info$iny
  ts <- nc_info$ts
  extraGlobalAttributes <- nc_info$extraGlobalAttributes
  if (fileType == "xy") proj4varname <- nc_info$proj4varname
  if (fileType == "xy") indataatt <- nc_info$indataatt
  if (fileType == "xy") ingridmapping <- nc_info$ingridmapping
  
  
  
  ## Remove the output file first (if there is one...)
  if (file.exists(outFile)) {
    invisible(file.remove(outFile))
  }
  
  outTimes <- seq(as.POSIXct(ts[1]), by = paste(nhourly, "hours"), length = (24 / nhourly) * length(ts))
  if (shiftouthours != 0) {
    outTimes <- outTimes + hours(shiftouthours)
  }
  outTatt <- paste0("hours since ", outTimes[1])
  outTvals <- seq(from = 0, length.out = length(outTimes), by = nhourly)
  
  if (fileType == "latlon") {
    dimX <- ncdim_def(name='lon', units='degrees_east', longname='longitude', vals=inx )
    dimY <- ncdim_def(name='lat', units='degrees_north', longname='latitude', vals=iny )
  } else if (fileType == "xy") {
    dimX <- ncdim_def(name='x', units='Meter', longname='x coordinate of projection', vals=inx )
    dimY <- ncdim_def(name='y', units='Meter', longname='y coordinate of projection', vals=iny )
  } else if (fileType == "curvilinear_2d") {
    print("add curv")
    dimX <- ncdim_def(name='x', units='', vals = 1:nx, create_dimvar=F)
    dimY <- ncdim_def(name='y', units='', vals = 1:ny, create_dimvar=F )
    dimNV4 <- ncdim_def(name='nv4', units='', vals = 1:4, create_dimvar=F )
  }
  dimT <- ncdim_def(name='time', units=outTatt, vals=outTvals )
  
  varData <- ncvar_def(name=outVar, units='-', dim=list(dimX,dimY,dimT), missval=NA, prec='double', compression = 5)
  
  if (fileType == "curvilinear_2d") {
    print("add curv")
    varLat <- ncvar_def(name='lat', units='degrees_north', dim=list(dimX,dimY), missval=NA, prec='double')
    varLon <- ncvar_def(name='lon', units='degrees_east', dim=list(dimX,dimY), missval=NA, prec='double')
    varLatBnds <- ncvar_def(name='lat_bnds', units='degrees_north', dim=list(dimX,dimY,dimNV4), missval=NA, prec='double')
    varLonBnds <- ncvar_def(name='lon_bnds', units='degrees_east', dim=list(dimX,dimY,dimNV4), missval=NA, prec='double')
    ncid_out <- nc_create( outFile, list(varLat,varLon,varLatBnds,varLonBnds,varData) )
  } else if (fileType == "xy") {
    varCoord <- ncvar_def(name=proj4varname, dim=list(), units="", prec='integer')
    ncid_out <- nc_create( outFile, list(varCoord, varData) )
  } else {
    ncid_out <- nc_create( outFile, list(varData) )
  }
  nc_close(ncid_out)
  
  ## Add data to file
  ncid_out <- nc_open( outFile, write = T )
  
  if (fileType == "curvilinear_2d") {
    print("add curv")
    ncvar_put(ncid_out, varLon, vals = inx)
    ncvar_put(ncid_out, varLat, vals = iny)
  } else if (fileType == "xy") {
    ## add PROJ4 attributes
    for (i in 1:length(ingridmapping$name)) {
      ncatt_put(ncid_out, proj4varname, ingridmapping$name[[i]], ingridmapping$value[[i]])
    }
  } else {
  }
  ncatt_put(ncid_out, "time", "calendar", "standard")
  
  ## Add variable attributes
  if (outVar == "tair") {
    for (i in 1:length(indataatt)) {
      ncatt_put(ncid_out, outVar, names(indataatt[i]), indataatt[[i]])
    }
    ncatt_put(ncid_out, outVar, "standard_name", "T")
    ncatt_put(ncid_out, outVar, "long_name", paste0(nhourly, " hourly temperature"))
  } else if (outVar == "potrad") {
    ncatt_put(ncid_out, outVar, "standard_name", "potrad")
    ncatt_put(ncid_out, outVar, "long_name", "potential radiation")
    ncatt_put(ncid_out, outVar, "units", "W m-2")
    if (fileType == "xy") ncatt_put(ncid_out, outVar, "grid_mapping", indataatt[["grid_mapping"]])
    if (fileType == "xy") ncatt_put(ncid_out, outVar, "esri_pe_string", indataatt[["esri_pe_string"]])
  } else {
    for (i in 1:length(indataatt)) {
      ncatt_put(ncid_out, outVar, names(indataatt[i]), indataatt[[i]])
    }
  }
  
  ## Add global attributes
  if (length(extraGlobalAttributes$name) > 0) {
    for (i in length(extraGlobalAttributes$name)) {
      ncatt_put(ncid_out, 0, extraGlobalAttributes$name, extraGlobalAttributes$att)
    }
  }
  
  nc_close(ncid_out)
}