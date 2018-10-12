#' @export
metGenRun <- function() {

  mgcheckVariables()
  
  nx <- metGen$settings$nx
  ny <- metGen$settings$ny

  ## makeOutputNetCDF
  makeNetcdfOut()
  
  ## DEFINE OUTPUT ARRAY
  outData <- NULL
  for (var in names(metGen$settings$outVars)) {
    outData[[var]] <- array(NA, dim = c(nx, ny, metGen$derived$nOutStepDay))
  }
  
  ### THE MAIN LOOP
  profile<-NULL
  profile$start.time.total <- Sys.time()
  for (iday in 1:metGen$derived$nday) {
    metGen$current$timestep <- (metGen$derived$nInStepDay * (iday-1)) + 1
    yday            <- metGen$derived$inYDays[iday]
    printf("Day: %d, date: %s\n", iday, metGen$derived$inDates[metGen$current$timestep])
    
    ## LOAD WHOLE DOMAIN FROM NETCDF
    profile$start.time.read <- Sys.time()
    inData <- readAllForcing(metGen$derived$inDates[metGen$current$timestep])
    profile$end.time.read <- Sys.time()
    
    # /*************************************************
    #   Precipitation
    # *************************************************/
    if(!is.null(metGen$settings$inVar$pr) && !is.null(outData$pr)) {
      nInStep <- metGen$derived$nInStepDay
      nOutStep <- metGen$derived$nOutStepDay
      outData$pr[]<-0
      if (nInStep <= nOutStep) {
        print("nOutStep is higher")
        inrec <- 1
        for(outrec in 1:nOutStep) {
          outData$pr[, , outrec] <- inData$pr[, , inrec]
          print(paste(outrec, inrec))
          if (!outrec%%(nOutStep/nInStep)) inrec <- inrec + 1
        }
      } else {
        print("nInStep is higher")
        outrec <- 1
        for(inrec in 1:nInStep) {
          outData$pr[, , outrec] <- outData$pr[, , outrec] + inData$pr[, , inrec]
          print(paste(outrec, inrec))
          if (!inrec%%(nInStep/nOutStep)) {
            outData$pr[, , outrec] <- outData$pr[, , outrec] / (nInStep/nOutStep)
            outrec <- outrec + 1
          }
        }
      }
    }
    # /*************************************************
    #   Shortwave radiation
    # *************************************************/
    radfrac <- rad_map_final_cr(metGen$derived$nOutStepDay, yday, gmt_float = 0, metGen$settings$lonlatbox)
    
    if(!is.null(metGen$settings$inVar$shortwave) && !is.null(outData$shortwave)) {
      for(rec in 1:metGen$derived$nOutStepDay) {
        outData$shortwave[, , rec] <- radfrac[ , , rec] * inData$shortwave[, , 1]
      }
      # print(image(outData$shortwave[, , 1]))
    }
    # /*************************************************
    #   Longwave radiation
    # *************************************************/
    if(!is.null(metGen$settings$inVar$longwave) && !is.null(outData$longwave)) {
      for(rec in 1:metGen$derived$nOutStepDay) {
        outData$longwave[, , rec] <- inData$longwave[, , 1]
      }
      
    }
    
    # /*************************************************
    #   Wind
    # *************************************************/
    if(!is.null(metGen$settings$inVar$wind) && !is.null(outData$wind)) {
      for(rec in 1:metGen$derived$nOutStepDay) {
        outData$wind[, , rec] <- inData$wind[, , 1]
      }
    }
    
    # /**************************************
    #   Atmospheric Pressure (Pa)
    # **************************************/
    if(!is.null(metGen$settings$inVar$pressure) && !is.null(outData$pressure)) {
      for(rec in 1:metGen$derived$nOutStepDay) {
        outData$pressure[, , rec] <- inData$pressure[, , 1]
      }
    }
    
    # /**************************************
    #   Temperature
    # **************************************/
    if(!is.null(metGen$settings$inVar$tasmin) && !is.null(metGen$settings$inVar$tasmin) && !is.null(outData$tas)) {
      outData$tas <- set_max_min_lonlat_cr(inData$tasmin[,,1], inData$tasmax[,,1], yday, metGen$derived$nOutStepDay, metGen$settings$lonlatbox)
    }
    
    # /*************************************************
    #   Vapor pressure
    # *************************************************/
    if(!is.null(metGen$settings$inVar$relhum) && !is.null(outData$vp)) {
      for(rec in 1:metGen$derived$nOutStepDay) {
        outData$vp <- set_vp_cr(outData$tas, inData$relhum[,,1], nx, ny, metGen$derived$nOutStepDay)
      }
    }
    if(!is.null(metGen$settings$inVar$qair) && !is.null(outData$vp)) {
      for(rec in 1:metGen$derived$nOutStepDay) {
        # outData$vp[,,rec] <- inData$qair[,,1] * inData$pressure[,,1]  / metGen$constants$EPS
        # outData$vp[,,rec] <- mg_sh2vp(inData$qair[,,1], outData$tas[,,rec], inData$pressure[,,1])
        # outData$vp[,,rec] <- sh2vp(inData$qair[,,1], inData$pressure[,,1])
        outData$vp[,,rec] <- sh2vp(inData$qair[,,1], inData$pressure[,,1])
      }
    }
    
    ## Convert to desired output unit
    for (var in names(metGen$settings$outVars)) {
        outData[[var]] <- convertUnit(outData[[var]], metGen$metadata$outvars[[var]]$internal_units, metGen$metadata$outvars[[var]]$output_units)
    }

    ## ADD OUTPUT TO NETCDF
    for (var in names(metGen$settings$outVars)) {
      timeIndex <- metGen$derived$nOutStepDay*(iday-1)+1
      metGen$settings$outVars[[var]]$ncid <- nc_open(metGen$settings$outVars[[var]]$filename, write = TRUE)
      ncvar_put(metGen$settings$outVars[[var]]$ncid,
                var,
                outData[[var]][,,],
                start = c(1, 1, timeIndex),
                count = c(nx, ny, metGen$derived$nOutStepDay)
      )
      nc_close(metGen$settings$outVars[[var]]$ncid)
    }
    rm(inData)
  }
  
  ncellsTotal <- nx*ny
  profile$end.time.total <- Sys.time()
  totTime<-as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")
  cat(sprintf("  Total: %.1f seconds for %d day(s). So: 100 years will take: %.2f days\n", 
              totTime, 
              metGen$derived$nday, 
              (totTime * 365.25 * 100) / metGen$derived$nday / 86400
  ))
}