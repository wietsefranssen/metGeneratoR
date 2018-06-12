mtclim_run<-function(settings = settings) {
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
  
  profile<-NULL
  
  ## makeOutputNetCDF
  for (iYear in 1:length(settings$ncOut)) {
    makeNetcdfOut(settings, elevation, settings$ncOut[[iYear]])
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
  el <- array(NA, dim = c(nx, ny, settings$intern$nrec_out))
  toNetCDFData <- list(el)[rep(1,length(settings$outputVars))]
  rm(el)
  
  ## LOAD WHOLE DOMAIN FROM NETCDF
  profile$start.time.read <- Sys.time()
  forcing_dataRTotal <- readAllForcing(settings, elevation)
  profile$end.time.read <- Sys.time()
  
  ## Init progressbar
  pb <- txtProgressBar(min = 0, max = ny, initial = 0, char = ">",
                       width = 80, title, label, style = 1, file = "")
  
  profile$start.time.run <- Sys.time()
  ## CELL LOOP
  for (iy in 1:ny) {
    output<-foreach(ix = 1:nx) %dopar% {
      if (!is.na(elevation$Data[ix,iy])) {
        settings$mtclim$elevation <- elevation$Data[ix,iy]
        settings$mtclim$lon<-elevation$xyCoords$x[ix]
        settings$mtclim$lat<-elevation$xyCoords$y[iy]
        
        ## RUN MLTCLIM
        # mtclimRun(forcing_dataR = selectForcingCell(settings, forcing_dataRTotal, ix, iy),
        #           settings = settings$mtclim)$out_data
        start<-1
        end<-length(settings$outputVars)*settings$intern$nrec_out
        c(start:end)
      }
    }
    
    for (ix in 1:length(elevation$xyCoords$x)) {
      if (!is.na(elevation$Data[ix,iy])) {
        ## ADD TO OUTPUT ARRAY
        for (iVar in 1:length(settings$outputVars)) {
          iStart <- ((iVar-1)*settings$intern$nrec_out)+1
          iEnd <- iVar*settings$intern$nrec_out
          toNetCDFData[[iVar]][ix,iy,] <- output[[ix]][iStart:iEnd]
        }
      }
    }
    rm(output)
    
    ## refresh progressbar
    setTxtProgressBar(pb, iy)
    
  }
  # Close ProgressBar
  close(pb)
  
  profile$end.time.run <- Sys.time()
  
  ## ADD OUTPUT TO NETCDF
  profile$start.time.write <- Sys.time()
  
  ## Define numer of years
  profile$start.time.write <- Sys.time()
  for (iYear in 1:length(settings$ncOut)) {
    ncid <- nc_open(settings$ncOut[[iYear]]$fileName, write = TRUE)
    sIndex <- settings$ncOut[[iYear]]$sIndex
    eIndex <- settings$ncOut[[iYear]]$eIndex
    for (iVar in 1:length(settings$outputVars))
    {
      ncvar_put(ncid,
                names(settings$outputVars)[iVar],
                toNetCDFData[[iVar]][,,sIndex:eIndex],
                start = c(1,
                          1,
                          1),
                count = c(nx,
                          ny,
                          settings$ncOut[[iYear]]$nrec_out)
      )
    }
    nc_close(ncid)
  }
  profile$end.time.write <- Sys.time()
  
  ## Print info about part
  cat(sprintf("  Times (read/run/write/total): %.1f/%.1f/%.1f/%.1f minutes\n",
              as.numeric(profile$end.time.read - profile$start.time.read, units = "mins"),
              as.numeric(profile$end.time.run - profile$start.time.run, units = "mins"),
              as.numeric(profile$end.time.write - profile$start.time.write, units = "mins"),
              as.numeric(profile$end.time.write - profile$start.time.read, units = "mins"),
              format(object.size(forcing_dataRTotal), units = "auto")))
  cat(sprintf("  Sizes (read/write): %s/%s\n",
              format(object.size(forcing_dataRTotal), units = "auto"),
              format(object.size(toNetCDFData), units = "auto")))
  
  
  cat(sprintf("\nFinished in %.1f minutes\n",as.numeric(Sys.time() - start.time.total, units = "mins")))
}
