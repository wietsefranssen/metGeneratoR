melo_get_nc_info <- function() {
  ncid_in <- nc_open(inFile, write = F)
  if (fileType == "latlon") {
    inx <- ncvar_get(ncid_in,"lon")
    nx <- length(inx)
    iny <- ncvar_get(ncid_in,"lat")
    ny <- length(iny)
    ## Make d2 lon and lat arrays
    ## First based on 1d lon lat files:
    lon_1d <- inx
    lat_1d <- iny
    lon_2d <- array(data = NA, dim = c(nx,ny))
    lat_2d <- array(data = NA, dim = c(nx,ny))
    for (ix in 1:nx) lon_2d[ix,] <- lon_1d[ix]  
    for (iy in 1:ny) lat_2d[,iy] <- lat_1d[iy]  
  } else if (fileType == "curvilinear_2d") {
    lon_2d <- inx <- ncvar_get(ncid_in,"lon")
    nx <- ncid_in$dim$x$len
    #inxatt <- ncatt_get(ncid_in,"lon")
    lat_2d <- iny <- ncvar_get(ncid_in,"lat")
    ny <- ncid_in$dim$y$len
    #inyatt <- ncatt_get(ncid_in,"lat")
  } else if (fileType == "xy") {
    inx <- ncvar_get(ncid_in,"x")
    nx <- ncid_in$dim$x$len
    iny <- ncvar_get(ncid_in,"y")
    ny <- ncid_in$dim$y$len
    indataatt <- ncatt_get(ncid_in,varname)
    ingridmapping <- nc_grid_mapping_atts(inFile)
    proj4string <- suppressWarnings(nc_gm_to_prj(ingridmapping))
    proj4varname <- ingridmapping$variable[1]
    x_2d <- array(data = NA, dim = c(nx,ny))
    y_2d <- array(data = NA, dim = c(nx,ny))
    for (ix in 1:nx) x_2d[ix,] <- inx[ix]  
    for (iy in 1:ny) y_2d[,iy] <- iny[iy]  
    x_1d <- as.vector(x_2d)
    y_1d <- as.vector(y_2d)
    xy <- data.frame(x=x_1d, y=y_1d)
    
    ## Transformed data
    pj <- proj4::project(xy, proj4string, inverse=TRUE)
    pj_lon <- pj$x
    dim(pj_lon) <- c(nx,ny)
    pj_lat <- pj$y
    dim(pj_lat) <- c(nx,ny)
    lat_2d <- pj_lat
    lon_2d <- pj_lon
  } else {
    stop("fileType not regonised")
  }
  intimes <- ncvar_get(ncid_in,"time")
  ts <- nc.get.time.series(ncid_in)
  if (shiftinhours != 0) {
    ts <- as.POSIXct(ts) + hours(shiftinhours)
  }
  
  # indata <- ncvar_get(ncid_in, varname)
  nc_close(ncid_in)

  extraGlobalAttributes <- NULL
  extraGlobalAttributes$name <- NULL
  extraGlobalAttributes$att <- NULL
  
  ## Checks
  if (shiftouthours != 0) {
    strmessage <- paste0("Output time shifted by ", shiftouthours, " hours")
    ## Print message
    cat(paste0("NOTE: ", strmessage, "\n"))
    ## Add global attribute
    ii <- length(extraGlobalAttributes$name) + 1
    extraGlobalAttributes$name[ii] <- paste0("disaggr_note_", ii)
    extraGlobalAttributes$att[ii] <- strmessage
  }
  if (shiftinhours != 0) {
    strmessage <- paste0("Input time shifted by ", shiftinhours, " hours")
    ## Print message
    cat(paste0("NOTE: ", strmessage, "\n"))
    ## Add global attribute
    ii <- length(extraGlobalAttributes$name) + 1
    extraGlobalAttributes$name[ii] <- paste0("disaggr_note_", ii)
    extraGlobalAttributes$att[ii] <- strmessage
  }
  
  ## Create function output
  nc_info <- NULL
  nc_info$nx <- nx
  nc_info$ny <- ny
  nc_info$lon_2d <- lon_2d
  nc_info$lat_2d <- lat_2d
  nc_info$inx <- inx
  nc_info$iny <- iny
  nc_info$ts <- ts
  nc_info$extraGlobalAttributes <- extraGlobalAttributes
  if (fileType == "xy") nc_info$proj4varname <- proj4varname
  if (fileType == "xy") nc_info$indataatt <- indataatt
  if (fileType == "xy") nc_info$ingridmapping <- ingridmapping
  
  return(nc_info)
}
