initSettings <- function(startdate = NULL,
                         enddate = NULL,
                         outstep = 24,
                         lonlatbox = NULL,
                         outfile = "output.nc",
                         outperyear = FALSE,
                         chunksizes = c(720,360,1)) {

  system <- list (
    nCores = 2,
    maxMem = 8  # in Gb
  );
  intern <- list (
    nrec_in = NULL,
    nrec_out = NULL
  );
  mtclim <- list (
    outstep = outstep,
    startyear = NULL,
    startmonth = NULL,
    startday = NULL,
    istartdate = NULL,
    endyear = NULL,
    endmonth = NULL,
    endday = NULL,
    ienddate = NULL,
    nForcing = NULL,
    forcingIds = NULL,
    alma_input = FALSE
  )
  
  settings <- list(
    lonlatbox = lonlatbox,
    startdate = as.Date(startdate),
    enddate = as.Date(enddate),
    outstep = NULL,
    outstep_per_day = NULL,
    outfile = outfile,
    outperyear = outperyear,
    chunksizes = chunksizes,
    mtclim = mtclim,
    intern = intern,
    system = system
  )
  
  settings <- setOutstep(outstep, settings)
  
  settings$intern$nrec_in = as.numeric((settings$enddate - settings$startdate) +1)
  settings$intern$nrec_out = settings$intern$nrec_in * (24/outstep)

  ## Add forcing info to settings
  nForcing = length(settings$inputVars)
  forcingIds = array(NA, dim = c(nForcing))
  for (iFor in 1:nForcing) {
    settings$mtclim$forcingIds[iFor] = settings$inputVars[[iFor]]$vicIndex
  }

  ymdStart<-as.numeric(strsplit(as.character(settings$startdate), "-")[[1]])
  ymdEnd<-as.numeric(strsplit(as.character(settings$enddate), "-")[[1]])

  settings$mtclim$startyear = ymdStart[1]
  settings$mtclim$startmonth = ymdStart[2]
  settings$mtclim$startday = ymdStart[3]
  settings$mtclim$istartdate = 1
  settings$mtclim$endyear = ymdEnd[1]
  settings$mtclim$endmonth = ymdEnd[2]
  settings$mtclim$endday = ymdEnd[3]
  settings$mtclim$ienddate = settings$intern$nrec_in

  ## make possible for splitup per year
  settings <- setNcOutInfo(settings)

  return(settings)
}

setOutstep <- function(outstep, settings = settings) {
  settings$outstep <- outstep
  settings$outstep_per_day <- 24 / outstep
  return(settings)
}

setInputVars <- function(settings, inputVars) {
  settings$inputVars <- inputVars

  if ( exists("alma", where = settings$inputVars[[1]]))  settings$mtclim$alma_input = settings$inputVars[[1]]$alma

  nForcing <- length(inputVars)
  forcingIds = array(NA, dim = c(nForcing))
  for (iFor in 1:nForcing) {
    forcingIds[iFor] = inputVars[[iFor]]$vicIndex
  }
  settings$mtclim$nForcing = nForcing
  settings$mtclim$forcingIds = forcingIds

  return(settings)
}

setNcOutInfo <- function(settings = settings) {
  settings$ncOut <- NULL
  ncOut <- NULL
  ncOut$startdate <- settings$startdate
  ncOut$enddate <- settings$enddate
  # ncOut$nrec_out <- settings$intern$nrec_out

  if (settings$outperyear) {
    ## Define numer of years
    sYear<-as.numeric(strsplit(as.character(settings$startdate), "-")[[1]])[1]
    eYear<-as.numeric(strsplit(as.character(settings$enddate), "-")[[1]])[1]
    nYear <- eYear - sYear + 1

    for (iYear in 1:nYear) {
      if (iYear == 1) {
        ncOut$startdate <- settings$startdate
        ncOut$sIndex <- 1
      } else {
        ncOut$startdate <- as.Date(paste0(sYear+iYear-1, "-01-01"))
        ncOut$sIndex <- settings$ncOut[[iYear - 1]]$sIndex + settings$ncOut[[iYear - 1]]$nrec_out
      }

      if (iYear == nYear) {
        ncOut$enddate <- settings$enddate
      } else {
        ncOut$enddate <- as.Date(paste0(sYear+iYear-1, "-12-31"))
      }
      nDays <-  as.numeric((ncOut$enddate - ncOut$startdate) + 1)
      ncOut$nrec_out <- nDays * (24/settings$outstep)
      ncOut$year <- as.numeric(strsplit(as.character(ncOut$startdate), "-")[[1]])[1]
      ncOut$eIndex <- (ncOut$sIndex + ncOut$nrec_out ) - 1
      ncOut$fileName <- paste0(settings$outfile, as.character(ncOut$year), ".nc")
      ## Fill list
      settings$ncOut[[iYear]] <- ncOut
    }
  } else {
    ncOut$startdate <- settings$startdate
    ncOut$sIndex <- 1
    ncOut$enddate <- settings$enddate
    nDays <-  as.numeric((ncOut$enddate - ncOut$startdate) + 1)
    ncOut$nrec_out <- nDays * (24/settings$outstep)
    ncOut$year <- as.numeric(strsplit(as.character(ncOut$startdate), "-")[[1]])[1]
    ncOut$eIndex <- (ncOut$sIndex + ncOut$nrec_out ) - 1
    ncOut$fileName <- paste0(settings$outfile,
                             gsub("-", "",as.character(ncOut$startdate)), "-",
                             gsub("-", "",as.character(ncOut$enddate)),
                             ".nc")
    settings$ncOut[[1]] <- ncOut
  }

  # print(settings$ncOut)
  return(settings)
}

makeNetcdfOut <- function(settings, mask, ncOut) {

  ## CREATE NETCDF
  FillValue <- 1e20

  ## Define dimensions
  dimX <- ncdim_def("lon", "degrees_east", mask$xyCoords$x)
  dimY <- ncdim_def("lat", "degrees_north",mask$xyCoords$y)
  timeString <-format(strptime(ncOut$startdate, format = "%Y-%m-%d", tz = "GMT"),format="%Y-%m-%d %T")
  timeArray <-c(0:(ncOut$nrec_out-1)) * (24 / (24/settings$outstep))
  dimT <- ncdim_def("time", paste0("hours since ",timeString), timeArray, unlim = FALSE)

  chunksizes_preffered <- settings$chunksizes
  dimsizes<-c(length(mask$xyCoords$x),length(mask$xyCoords$y),ncOut$nrec_out)
  chunksizes<- pmin(chunksizes_preffered,dimsizes)
  ################
  # data <- ncvar_def(name=names(settings$outputVars[1]), units='', dim=list(dimX,dimY,dimT), missval=FillValue, prec="float")
  data <- ncvar_def(name=names(settings$outputVars[1]), units='', compression = 7, chunksizes=chunksizes, dim=list(dimX,dimY,dimT), missval=FillValue, prec="float")
  dataAllVars <- list(data)[rep(1,length(settings$outputVars))]
  for (iVar in 1:length(settings$outputVars))
  {
    dataAllVars[[iVar]] <- ncvar_def(name=names(settings$outputVars[iVar]), units='', compression = 7, chunksizes=chunksizes, dim=list(dimX,dimY,dimT), missval=FillValue, prec="float")
    # dataAllVars[[iVar]] <- ncvar_def(name=names(settings$outputVars[iVar]), units='', dim=list(dimX,dimY,dimT), missval=FillValue, prec="float")
  }

  ## SAVE AS NC-DATA
  cat(sprintf("Output files: %s\n", ncOut$fileName))
  # ncid <- nc_create_par(settings$outfile, dataAllVars,force_v4=TRUE)
  ncid <- nc_create(ncOut$fileName, dataAllVars, force_v4=TRUE)

  ncatt_put( ncid, "lon", "standard_name", "longitude")
  ncatt_put( ncid, "lon", "long_name",     "Longitude")
  ncatt_put( ncid, "lon", "axis",          "X")
  ncatt_put( ncid, "lat", "standard_name", "latitude")
  ncatt_put( ncid, "lat", "long_name",     "Latitude")
  ncatt_put( ncid, "lat", "axis",          "Y")
  ncatt_put( ncid, "time", "standard_name", "time")
  ncatt_put( ncid, "time", "calendar",     "standard")
  for (iVar in 1:length(settings$outputVars))
  {
    ncatt_put( ncid, names(settings$outputVars[iVar]), "standard_name", names(settings$outputVars[iVar]))
    ncatt_put( ncid, names(settings$outputVars[iVar]), "long_name", settings$outputVars[[iVar]]$longName)
    ncatt_put( ncid, names(settings$outputVars[iVar]), "vic_name", settings$outputVars[[iVar]]$VICName)
    ncatt_put( ncid, names(settings$outputVars[iVar]), "units", settings$outputVars[[iVar]]$units)
  }
  # ## Global Attributes
  ncatt_put( ncid, 0, "NetcdfCreatationDate", as.character(Sys.Date()))

  ## Close Netcdf file
  nc_close(ncid)
}

calcMinNParts <- function(settings, mask) {
  BYTE<-8

  maxMem <- settings$system$maxMem * 1024^3

  ## Calculate memory needed:
  memInput <- length(mask$Data) * settings$intern$nrec_in * length(settings$inputVars) * BYTE
  memOutput <- length(mask$Data) * settings$intern$nrec_out * length(settings$outputVars) * BYTE
  memExtra <- length(mask$Data) * 100 * BYTE
  memTotal <- (memInput + memOutput + memExtra)

  minNParts <- ceiling(memTotal / maxMem) #
  # maxNParts <- ceiling(memTotal / maxMem * 2) #*2 for safety

  cat(sprintf("Memory needed: %s (Max mem: %s), minimum number of parts: %d\n", hsize(memTotal), hsize(maxMem), minNParts))

  return(minNParts)
}

setSubDomains <- function(settings, mask, nPart = NULL) {
  # nPart<-1
  # mask<-elevation

  nCells<-length(mask$xyCoords$x)*length(mask$xyCoords$y)
  nActive <- length(mask$Data[!is.na(mask$Data)])

  if (nPart == 1) {
    part <- list(sx = 1,
                 nx = length(mask$xyCoords$x),
                 sy = 1,
                 ny = length(mask$xyCoords$y),
                 slon = min(mask$xyCoords$x),
                 elon = max(mask$xyCoords$x),
                 slat = min(mask$xyCoords$y),
                 elat = max(mask$xyCoords$y))
    parts <- list(part)[rep(1,nPart)]
  } else {

    nxxx<-length(mask$xyCoords$x)
    nyyy<-length(mask$xyCoords$y)
    lengthxy<-ceiling(sqrt(nCells/nPart))
    parts<-NULL
    ey<-ex<-0
    sy<-sx<-0
    ny<-nx<-0
    sxPrev<-syPrev<-1
    i<-1
    while(ey < nyyy) {
      while(ex < nxxx) {
        sx[[i]]<-sxPrev
        sy[[i]]<-syPrev
        nx[[i]]<-lengthxy
        ny[[i]]<-lengthxy

        ex<-sxPrev+lengthxy-1
        sxPrev<-sx[[i]]+lengthxy
        i<-i+1
      }
      nx[[i-1]]<-nxxx-sx[[i-1]]+1
      ey<-syPrev+lengthxy-1
      syPrev<-syPrev+lengthxy
      ex<-1
      sxPrev<-1
    }
    nPart<-length(sy)

    for (iPart in 1:nPart) {
      if ((nyyy-(sy[[iPart]])) < lengthxy) {
        ny[[iPart]]<-nyyy- sy[[iPart]]+1
      }
      if ((nxxx-(sx[[iPart]])) < lengthxy) {
        nx[[iPart]]<-nxxx- sx[[iPart]]+1
      }
    }
    part <- list(sx = 1,
                 nx = 1,
                 sy = NULL,
                 ny = NULL,
                 slon = min(mask$xyCoords$x),
                 elon = max(mask$xyCoords$x),
                 slat = NULL,
                 elat = NULL)
    parts <- list(part)[rep(1,nPart)]

    ## The devision
    for (iPart in 1:nPart) {
      endx<- sx[iPart]+nx[iPart]-1
      endy<- sy[iPart]+ny[iPart]-1
    }

    for (iPart in 1:(nPart)) {
      endx<- sx[iPart]+nx[iPart]-1
      endy<- sy[iPart]+ny[iPart]-1
      parts[[iPart]]$sx <- sx[iPart]
      parts[[iPart]]$sy <- sy[iPart]
      parts[[iPart]]$slon <- mask$xyCoords$x[parts[[iPart]]$sx]
      parts[[iPart]]$elon <- mask$xyCoords$x[endx]
      parts[[iPart]]$slat <- mask$xyCoords$y[parts[[iPart]]$sy]
      parts[[iPart]]$elat <- mask$xyCoords$y[endy]
      parts[[iPart]]$nx <- nx[iPart]
      parts[[iPart]]$ny <- ny[iPart]
    }
  }

  cat(sprintf("Total cells: %d, Active Cells: %d, nParts: %d\n",nCells , nActive, nPart))

  return(parts)
}

readForcing <- function(settings, iPart) {
  forcing_dataR <- list()
  for (i in 1:length(settings$inputVars)) {
    forcing_dataR[[i]]<- array(0, dim=c(settings$intern$nrec_in,settings$parts[[iPart]]$ny,settings$parts[[iPart]]$nx))
  }
  ## Read data
  sx <- settings$parts[[iPart]]$sx
  sy <- settings$parts[[iPart]]$sy
  nx <- settings$parts[[iPart]]$nx
  ny <- settings$parts[[iPart]]$ny
  for (iVar in 1:length(settings$inputVars)) {
    forcing_dataR[[iVar]][] <- ncLoad(file = settings$inputVars[[iVar]]$ncFileName,
                                      varName = settings$inputVars[[iVar]]$ncName,
                                      lonlatbox = c(settings$parts[[iPart]]$slon,
                                                    settings$parts[[iPart]]$elon,
                                                    settings$parts[[iPart]]$slat,
                                                    settings$parts[[iPart]]$elat),
                                      timesteps = c(1:settings$intern$nrec_in))$Data
  }
  return(forcing_dataR)
}

readForcingAll <- function(part, settings, mask) {
  # mask<-elevation
  forcing_dataR <- list()
  for (i in 1:length(settings$inputVars)) {
    forcing_dataR[[i]]<- array(0, dim=c(settings$intern$nrec_in,  part$ny,  part$nx))
  }
  ## Read data
  for (iVar in 1:length(settings$inputVars)) {
    forcing_dataR[[iVar]][] <- ncLoad(file = settings$inputVars[[iVar]]$ncFileName,
                                      varName = settings$inputVars[[iVar]]$ncName,
                                      # lonlatbox = c(settings$lonlatbox[1],
                                      #               settings$lonlatbox[2],
                                      #               settings$lonlatbox[3],
                                      #               settings$lonlatbox[4]),
                                      lonlatbox = c(part$slon,
                                                    part$elon,
                                                    part$slat,
                                                    part$elat),
                                      timesteps = c(1:settings$intern$nrec_in))$Data
  }
  return(forcing_dataR)
}

readAllForcing <- function(settings, mask, timestep) {
  # mask<-elevation
  nx <- length(mask$xyCoords$x)
  ny <- length(mask$xyCoords$y)
  
  forcing_dataR <- list()
  for (i in 1:length(settings$inputVars)) {
    forcing_dataR[[i]]<- array(0, dim=c(nx, ny, 1))
  }
  ## Read data
  for (iVar in 1:length(settings$inputVars)) {
    forcing_dataR[[iVar]][] <- ncLoad(file = settings$inputVars[[iVar]]$ncFileName,
                                      varName = settings$inputVars[[iVar]]$ncName,
                                      # lonlatbox = c(settings$lonlatbox[1],
                                      #               settings$lonlatbox[2],
                                      #               settings$lonlatbox[3],
                                      #               settings$lonlatbox[4]),
                                      lonlatbox = settings$lonlatbox,
                                      timesteps = timestep)$Data
  }
  return(forcing_dataR)
}

selectForcingCell <-function(settings, forcing_dataRTotal, ix, iy) {
  forcing_dataR<-NULL

  for (i in 1:24) {
    forcing_dataR[[i]]<- array(0, dim=c(settings$intern$nrec_in))
  }
  ## Read data
  for (iVar in 1:length(settings$inputVars)) {
    forcing_dataR[[settings$inputVars[[iVar]]$vicIndex+1]][] <- forcing_dataRTotal[[iVar]][,iy,ix]
  }
  return(forcing_dataR)
}

initParams <- function(startdate = NULL, enddate = NULL, outstep = 24, lonlatbox = NULL) {
  intern <- list (
    lonlatboxindex = NULL,
    lons = NULL,
    lats = NULL,
    nx = NULL,
    ny = NULL,
    startx = NULL,
    starty = NULL,
    nrec_in = NULL,
    nrec_out = NULL
  );
  toCpp <- list (
    outstep = outstep,
    startyear = NULL,
    startmonth = NULL,
    startday = NULL,
    istartdate = NULL,
    endyear = NULL,
    endmonth = NULL,
    endday = NULL,
    ienddate = NULL,
    nForcing = NULL,
    forcingIds = NULL,
    alma_input = FALSE
  )
  params <- list(
    lonlatbox = lonlatbox,
    startdate = as.Date(startdate),
    enddate = as.Date(enddate),
    outstep = outstep,
    toCpp = toCpp
  )
  return(params)
}

readForcingNetCDF <- function(variableInfo, params) {
  ## Read data
  nx <- params$intern$nx
  ny <- params$intern$ny
  startx <- params$intern$startx
  starty <- params$intern$starty
  forcing_dataR <- list()
  for (i in 1:24) {
    forcing_dataR[[i]]<- array(0, dim=c(nx,ny,params$intern$nrec_in))
  }
  for (iVar in 1:length(variableInfo)) {
    ncid <- nc_open(variableInfo[[iVar]]$ncFileName)
    forcing_dataR[[variableInfo[[iVar]]$vicIndex+1]][] <- ncvar_get(nc = ncid,
                                                                    varid =  variableInfo[[iVar]]$ncName,
                                                                    start = c(startx,
                                                                              starty,
                                                                              params$toCpp$istartdate),
                                                                    count = c(nx,ny,params$intern$nrec_in)
    );
    nc_close(ncid)

    ## Reverse the lat array
    if (params$intern$lonlatboxindex[4] > params$intern$lonlatboxindex[3]) {
      # forcing_dataR[[variableInfo[[iVar]]$vicIndex+1]][] <- apply(forcing_dataR[[variableInfo[[iVar]]$vicIndex+1]][], 3, rev)
    }
  }


  return(forcing_dataR)
}

readElevationNetCDF <- function(elevation, params) {
  ## Read data
  nx <- params$intern$nx
  ny <- params$intern$ny
  startx <- params$intern$startx
  starty <- params$intern$starty
  elevation$data <- array(0, dim=c(nx,ny))
  ncid <- nc_open(elevation$ncFileName)
  elevation$data[] <- ncvar_get(nc = ncid,
                                varid =  elevation$ncName,
                                start = c(startx,
                                          starty),
                                count = c(nx,ny)
  );
  nc_close(ncid)


  return(elevation)
}
