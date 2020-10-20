#' @export
metGenRun <- function() {
  
  mgcheckVariables()
  
  ## makeOutputNetCDF
  makeNetcdfOut()
  
  ## DEFINE OUTPUT ARRAY
  outData <- NULL
  for (var in names(metGen$settings$outVars)) {
    outData[[var]] <- array(NA, dim = c(metGen$settings$nx, metGen$settings$ny, metGen$derived$nOutStepDay))
  }
  
  nInStep <- metGen$derived$nInStepDay
  nOutStep <- metGen$derived$nOutStepDay
  
  maxStep <- max(nInStep, nOutStep)
  if (nOutStep > nInStep) {
    inrecs <- rep(1:nInStep, each = ceiling(nOutStep/nInStep))
    outrecs <- rep(1:nOutStep)
  } else {
    inrecs <- rep(1:nInStep)
    outrecs <- rep(1:nOutStep, each = ceiling(nInStep/nOutStep))
  }
  
  ### THE MAIN LOOP
  profile<-NULL
  profile$start.time.total <- Sys.time()
  for (iday in 1:metGen$derived$nday) {
    # iday<-1
    metGen$current$timestep <- (metGen$derived$nInStepDay * (iday-1)) + 1
    yday            <- metGen$derived$inYDays[iday]
    printf("Day: %d, date: %s\n", iday, metGen$derived$inDates[metGen$current$timestep])
    
    ## LOAD WHOLE DOMAIN FROM NETCDF
    profile$start.time.read <- Sys.time()
    inData <- readAllForcing(metGen$derived$inDates[metGen$current$timestep])
    profile$end.time.read <- Sys.time()
    
    
    # /*************************************************
    #   radiation fraction
    # *************************************************/
    if (metGen$output$ndim == 1) {
      lonlat2d <- F 
    } else {
      lonlat2d <- T
    }
    if (!is.null(outData$radfrac)) {
      radfrac <- rad_map_final_2dll_cr(metGen$derived$nOutStepDay, yday, gmt_float = 0, 
                                              metGen$settings$xybox, 
                                              metGen$output$lats,
                                              lonlat2d)
      print(radfrac[1,1 , ])
      for(i in 1:maxStep) outData$radfrac[, ,outrecs[i]] <- radfrac[, , outrecs[i]]
    }    
    
    # /*************************************************
    #   Precipitation
    # *************************************************/
    if(!is.null(outData$pr)) {
      if(metGen$metadata$inVars$pr$enabled) { ## if precip is provided
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) outData$pr[, ,outrecs[i]] <- inData$pr[, , inrecs[i]]
          # for(i in 1:maxStep) outData$pr[, ,outrecs[i]] <- inData$pr[, , 1]
        } else { ## aggregate to lower number of timesteps
          outData$pr[] <- 0
          for(i in 1:maxStep)  outData$pr[, , outrecs[i]] <- outData$pr[, , outrecs[i]] + ( inData$pr[, , inrecs[i]] / (nInStep/nOutStep) )
        }
      } else { ## if rainf and snowf are provided
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) outData$pr[, ,outrecs[i]] <- inData$rainf[, , inrecs[i]] + inData$snowf[, , inrecs[i]]
        } else { ## aggregate to lower number of timesteps
          outData$pr[] <- 0
          for(i in 1:maxStep) outData$pr[, , outrecs[i]] <- outData$pr[, , outrecs[i]] + ( (inData$rainf[, , inrecs[i]] + inData$snowf[, , inrecs[i]])  / (nInStep/nOutStep) )
        }
      }
    }
    
    # /*************************************************
    #   Shortwave radiation
    # *************************************************/
    if (!is.null(outData$swdown)) {
      if(metGen$metadata$inVars$swdown$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          # radfrac <- rad_map_final_2dll_cr(metGen$derived$nOutStepDay, yday, gmt_float = 0,
                                           # metGen$settings$xybox,
                                           # metGen$output$lats,
                                           # lonlat2d)
          # print(dim(radfrac))
          # print(dim(outData$swdown))
          # for(i in 1:maxStep) outData$swdown[, , outrecs[i]] <- radfrac[ , , outrecs[i]]
          lats <- aperm(array(metGen$output$lats, dim=c(360,720)),c(2,1))
          lons <- aperm(array(metGen$output$lons, dim=c(720,360)),c(1,2))
          radfrac <- rad_map_final_cr(metGen$derived$nOutStepDay, yday, gmt_float = 0, metGen$settings$xybox, lats, lons)
          print(dim(radfrac))
          print(dim(outData$swdown))
          image(radfrac[,,1])
          swdown_day <- apply(inData$swdown, c(1,2), mean)
          for(i in 1:maxStep) outData$swdown[, , outrecs[i]] <- radfrac[ , , outrecs[i]] #* swdown_day
        } else { ## aggregate to lower number of timesteps
          outData$swdown[]<-0
          for(i in 1:maxStep) outData$swdown[, , outrecs[i]] <- outData$swdown[, , outrecs[i]] + ( inData$swdown[, , inrecs[i]] /  (nInStep/nOutStep) )
        }
      }
    }    
    
    # /*************************************************
    #   Longwave radiation
    # *************************************************/
    if (!is.null(outData$lwdown)) {
      if(metGen$metadata$inVars$lwdown$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) outData$lwdown[, ,outrecs[i]] <- inData$lwdown[, , inrecs[i]]
        } else { ## aggregate to lower number of timesteps
          outData$lwdown[] <- 0
          for(i in 1:maxStep) outData$lwdown[, , outrecs[i]] <- outData$lwdown[, , outrecs[i]] + ( inData$lwdown[, , inrecs[i]] / (nInStep/nOutStep) )
        }
      }
    }
    
    # /*************************************************
    #   Wind
    # *************************************************/
    if (!is.null(outData$wind)) {
      if(metGen$metadata$inVars$wind$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) outData$wind[, ,outrecs[i]] <- inData$wind[, , inrecs[i]]
        } else { ## aggregate to lower number of timesteps
          outData$wind[] <- 0
          for(i in 1:maxStep) outData$wind[, , outrecs[i]] <- outData$wind[, , outrecs[i]] + ( inData$wind[, , inrecs[i]] / (nInStep/nOutStep) )
        }
      }
    }
    
    # /**************************************
    #   Atmospheric Pressure (Pa)
    # **************************************/
    if (!is.null(outData$psurf)) {
      if(metGen$metadata$inVars$psurf$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) outData$psurf[, ,outrecs[i]] <- inData$psurf[, , inrecs[i]]
        } else { ## aggregate to lower number of timesteps
          outData$psurf[] <- 0
          for(i in 1:maxStep) outData$psurf[, , outrecs[i]] <- outData$psurf[, , outrecs[i]] + ( inData$psurf[, , inrecs[i]] / (nInStep/nOutStep) )
        }
      }
    }
    
    # /**************************************
    #   Temperature
    # **************************************/
    if (!is.null(outData$tas)) {
      if(metGen$metadata$inVars$tasmin$enabled && metGen$metadata$inVars$tasmax$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          if (nInStep > 1) stop(printf("Dissagregation of \"tasmin\" and \"tasmax\" into \"tas\" is only possible for daily input!"))
          outData$tas <- set_max_min_lonlat_cr(inData$tasmin[,,1], inData$tasmax[,,1], yday, metGen$derived$nOutStepDay, metGen$settings$xybox)
        } else { ## aggregate to lower number of timesteps
          outData$tas[] <- 0
          for(i in 1:maxStep) outData$tas[, , outrecs[i]] <- outData$tas[, , outrecs[i]] + 
              ( ( inData$tasmin[, , inrecs[i]] + inData$tasmax[, , inrecs[i]] ) / ((nInStep/nOutStep)*2) )
        }
      } else { ## if tas in provided as input
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) outData$tas[, ,outrecs[i]] <- inData$tas[, , inrecs[i]]
        } else { ## aggregate to lower number of timesteps
          outData$tas[] <- 0
          for(i in 1:maxStep) outData$tas[, , outrecs[i]] <- outData$tas[, , outrecs[i]] + ( inData$tas[, , inrecs[i]] / (nInStep/nOutStep) )
        }
      }
    }
    
    # /**************************************
    #   Min and Max Temperature
    # **************************************/
    if (!is.null(outData$tasmin)) {
      if(metGen$metadata$inVars$tas$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          stop(printf("disaggregation of \"tas\" into \"tasmin\" is not possible!"))
        } else { ## aggregate to lower number of timesteps
          outData$tasmin[] <- 0
          for(i in 1:maxStep) outData$tasmin[, , outrecs[i]] <- pmin(outData$tasmin[, , outrecs[i]], inData$tas[, , inrecs[i]])
        }
      } else stop(printf("tas should be provided!\n"))
    }
    
    if (!is.null(outData$tasmax)) {
      if(metGen$metadata$inVars$tas$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          stop(printf("disaggregation of \"tas\" into \"tasmax\" is not possible!"))
        } else { ## aggregate to lower number of timesteps
          outData$tasmax[] <- 0
          for(i in 1:maxStep) outData$tasmax[, , outrecs[i]] <- pmax(outData$tasmax[, , outrecs[i]], inData$tas[, , inrecs[i]])
        }
      } else stop(printf("tas should be provided!\n"))
    }
    
    # /*************************************************
    #   Vapor pressure
    # *************************************************/
    # if(!is.null(metGen$settings$inVars$relhum) && !is.null(outData$vp)) {
    #   for(rec in 1:metGen$derived$nOutStepDay) {
    #     outData$vp <- set_vp_cr(outData$tas, inData$relhum[,,1], metGen$settings$nx, metGen$settings$ny, metGen$derived$nOutStepDay)
    #   }
    # }
    # if(!is.null(metGen$settings$inVars$qair) && !is.null(outData$vp)) {
    #   for(rec in 1:metGen$derived$nOutStepDay) {
    #     # outData$vp[,,rec] <- inData$qair[,,1] * inData$psurf[,,1]  / metGen$constants$EPS
    #     # outData$vp[,,rec] <- mg_sh2vp(inData$qair[,,1], outData$tas[,,rec], inData$psurf[,,1])
    #     # outData$vp[,,rec] <- sh2vp(inData$qair[,,1], inData$psurf[,,1])
    #     outData$vp[,,rec] <- sh2vp(inData$qair[,,1], inData$psurf[,,1])
    #   }
    # }
    if (!is.null(outData$vp)) {
      if (metGen$metadata$inVars$relhum$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) {
            outData$vp[, , ] <- set_vp_cr(outData$tas[, , ], inData$relhum[, , inrecs[i]], metGen$settings$nx, metGen$settings$ny, metGen$derived$nOutStepDay)
          }
        } else { ## aggregate to lower number of timesteps
          outData$vp[] <- 0
          for(i in 1:maxStep) outData$vp[, , outrecs[i]] <- outData$vp[, , outrecs[i]] + 
              ( set_vp_cr(outData$tas[, ,outrecs[i]], inData$relhum[, , inrecs[i]], metGen$settings$nx, metGen$settings$ny, metGen$derived$nOutStepDay) / 
                  (nInStep/nOutStep) )
        }
      } else if(metGen$metadata$inVars$qair$enabled && metGen$metadata$inVars$psurf$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) outData$vp[,,outrecs[i]] <- sh2vp(inData$qair[,,inrecs[i]], inData$psurf[,,inrecs[i]])
        } else { ## aggregate to lower number of timesteps
          outData$vp[] <- 0
          for(i in 1:maxStep) outData$vp[, , outrecs[i]] <- outData$vp[, , outrecs[i]] +  ( sh2vp(inData$qair[,,inrecs[i]], inData$psurf[,,inrecs[i]]) /  (nInStep/nOutStep) )
        }
      } else if(metGen$metadata$inVars$vp$enabled) {
        if (nInStep < nOutStep) { ## disaggregate to higher number of timesteps
          for(i in 1:maxStep) outData$vp[, ,outrecs[i]] <- inData$vp[, , inrecs[i]]
        } else { ## aggregate to lower number of timesteps
          outData$vp[] <- 0
          for(i in 1:maxStep) outData$vp[, , outrecs[i]] <- outData$vp[, , outrecs[i]] + ( inData$vp[, , inrecs[i]] / (nInStep/nOutStep) )
        }
      } else {
        if (nInStep > 1) stop(printf("Cannot convert to vp!"))
      }
    }
    
    ## Convert to desired output unit
    for (var in names(metGen$settings$outVars)) {
      outData[[var]] <- convertUnit(outData[[var]], 
                                    metGen$metadata$outVars[[var]]$internal_units, 
                                    metGen$metadata$outVars[[var]]$output_units,
                                    metGen$settings$outDt)
    }
    
    ## ADD OUTPUT TO NETCDF
    # need to switch dimension order first (fix later)
    # outData[[var]]<-aperm(outData[[var]],c(2,1,3))
    metGen$outData <- outData
    metGen$inData <- inData
    
    for (var in names(metGen$settings$outVars)) {
      timeIndex <- metGen$derived$nOutStepDay*(iday-1)+1
      # metGen$settings$outVars[[var]]$ncid <- nc_open(metGen$settings$outVars[[var]]$filename, write = TRUE)
      metGen$settings$outVars[[var]]$ncid <- open.nc(metGen$settings$outVars[[var]]$filename, write = TRUE)
      
      # ncvar_put(metGen$settings$outVars[[var]]$ncid,
      #           var,
      #           outData[[var]][,,],
      #           start = c(1, 1, timeIndex),
      #           count = c(metGen$settings$ny, metGen$settings$nx, metGen$derived$nOutStepDay)
      # )
      var.put.nc(metGen$settings$outVars[[var]]$ncid,
                 var,
                 outData[[var]][,,],
                 start = c(1, 1, timeIndex),
                 count = c(metGen$settings$nx, metGen$settings$ny, metGen$derived$nOutStepDay)
      )
      # nc_close(metGen$settings$outVars[[var]]$ncid)
      close.nc(metGen$settings$outVars[[var]]$ncid)
    }
    rm(inData)
  }
  
  ncellsTotal <- metGen$settings$nx*metGen$settings$ny
  profile$end.time.total <- Sys.time()
  totTime<-as.numeric(profile$end.time.total   - profile$start.time.total, units = "secs")
  cat(sprintf("  Total: %.1f seconds for %d day(s). So: 100 years will take: %.2f days\n", 
              totTime, 
              metGen$derived$nday, 
              (totTime * 365.25 * 100) / metGen$derived$nday / 86400
  ))
}
