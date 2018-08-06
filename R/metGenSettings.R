metGen <- new.env()

mgsetLonlatbox <- function(lonlatbox) {
  metGen$settings$lonlatbox <- lonlatbox
  metGen$settings$x <- seq(lonlatbox[1],lonlatbox[2], 0.5)
  metGen$settings$y <- seq(lonlatbox[3],lonlatbox[4], 0.5)
  metGen$settings$nx <- length(metGen$settings$x)
  metGen$settings$ny <- length(metGen$settings$y)
}

mgsetInDt <- function(inDt) {
  metGen$derived$inDt <- metGen$settings$inDt <- inDt
  metGen$derived$nInStepDay  <- 24 / metGen$derived$inDt
  
  ## Update mgsetPeriod
  if (!is.null(metGen$settings$startDate) && !is.null(metGen$settings$endDate))
  {
    mgsetPeriod(metGen$settings$startDate, metGen$settings$endDate)
  }
}

mgsetOutDt <- function(outDt) {
  metGen$derived$outDt <- metGen$settings$outDt <- outDt
  metGen$derived$nOutStepDay  <- 24 / metGen$derived$outDt
  
  ## Update mgsetPeriod
  if (!is.null(metGen$settings$startDate) && !is.null(metGen$settings$endDate))
  {
    mgsetPeriod(metGen$settings$startDate, metGen$settings$endDate)
  }
}

mgsetPeriod <- function(startdate, enddate) {
  metGen$settings$startDate <- startdate
  metGen$settings$endDate <- enddate
  metGen$derived$startDate <- as.POSIXlt(startdate, tz = "GMT")
  metGen$derived$endDate <- as.POSIXlt(enddate, tz = "GMT")
  
  
  metGen$derived$nday <- as.numeric(metGen$derived$endDate - metGen$derived$startDate + 1)
  metGen$derived$nrec_in <- metGen$derived$nday * metGen$derived$nInStepDay
  metGen$derived$nrec_out <- metGen$derived$nday * metGen$derived$nOutStepDay
  
  by <- paste(metGen$derived$inDt, "hours")
  metGen$derived$inDates = seq(metGen$derived$startDate, length = metGen$derived$nrec_in, by = by)
  metGen$derived$inYDays = as.POSIXlt(metGen$derived$inDates)$yday + 1
  by <- paste(metGen$derived$outDt, "hours")
  metGen$derived$outDates = seq(metGen$derived$startDate, length = metGen$derived$nrec_out, by = by)
  metGen$derived$outYDays = as.POSIXlt(metGen$derived$outDates)$yday + 1
}

mgsetElevation <- function(ncname, filename) {
  metGen$settings$elevation$ncname <- ncname
  metGen$settings$elevation$filename <- filename
}

mgsetInVars <- function(varlist) {
  metGen$settings$inVar <- NULL
  for (var in names(varlist)) {
    ## Check if input var exists
    if (!is.element(var, names(metGen$metadata$invars))) {
      cat(paste0("WARING: The variable ", var, " is not a standard input variable!\n",
                 "        Valid variables are: ", paste(names(metGen$metadata$invars),collapse=", ") , "\n",
                 "        Ingoring: ", var, "\n"))
    } else {
      metGen$settings$inVar[[var]]$ncname <- varlist[[var]]$ncname
      metGen$settings$inVar[[var]]$filename <- varlist[[var]]$filename
      metGen$settings$inVar[[var]]$longName <- metGen$metadata$invars[[var]]$longName
      # metGen$settings$inVar[[var]]$units <- metGen$metadata$invars[[var]]$units
    }
  }
  mgcheckInVars()
}

mgsetOutVars <- function(varnames) {
  ## Check if inputvars are already defined
  if (length(metGen$settings$inVar) <= 0) {
    cat(paste0("First set the input variables. Based on that the possible output variables can be defined\n"))
    stop()
  } else {
    metGen$settings$outVars <- NULL
    for (var in varnames) {
      ## Check input variables
      if (var == "shortwave" && is.null(metGen$settings$inVar[["shortwave"]])) {
        stop(paste0("Shortwave output can only be generated if shortwave is provided as input\nOn the mmoment, only the following input variables are defined: ", paste(names(metGen$settings$inVar), collapse=", "), "\n"), call. = FALSE)
      }
      if (var == "pr" && is.null(metGen$settings$inVar[["pr"]])) {
        stop(paste0("Precipitation output can only be generated if precipitation is provided as input\nOn the mmoment, only the following input variables are defined: ", paste(names(metGen$settings$inVar), collapse=", "), "\n"), call. = FALSE)
      }
      
      ## Check if output var exists
      if (!is.element(var, names(metGen$metadata$outvars))) {
        cat(paste0("WARING: The variable ", var, " is not a standard output variable!\n",
                   "        Valid variables are: ", paste(names(metGen$metadata$outvars),collapse=", ") , "\n",
                   "        Ingoring: ", var, "\n"))
      } else {
        metGen$settings$outVars[[var]]$name <- var
        metGen$settings$outVars[[var]]$longName <- metGen$metadata$outvars[[var]]$longName
        metGen$settings$outVars[[var]]$units <- metGen$metadata$outvars[[var]]$units
      }
    }
  }
  mgsetOutName(nameString = "output/<VAR>/<VAR>_<SYEAR><SMONTH><SDAY>_<EYEAR><EMONTH><EDAY>.nc", message = F)
}

mgsetOutName <- function(nameString, message = F) {
  # nameString <- "output/<VAR>/<VAR>_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_<SYEAR>.nc"
  
  for (var in names(metGen$settings$outVars)) {
    nameStringTmp <- nameString
    ## Substitute Dates
    nameStringTmp <- gsub("<SYEAR>", format(metGen$derived$startDate,format="%Y"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<EYEAR>", format(metGen$derived$endDate,format="%Y"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<SMONTH>", format(metGen$derived$startDate,format="%m"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<EMONTH>", format(metGen$derived$endDate,format="%m"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<SDAY>", format(metGen$derived$startDate,format="%m"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<EDAY>", format(metGen$derived$endDate,format="%m"), nameStringTmp, ignore.case = T)
    ## Substitute Variable name
    nameStringTmp <- gsub("<VAR>", var, nameStringTmp, ignore.case = T)
    
    metGen$settings$outVars[[var]]$filename <- nameStringTmp
  }
  
  ## print results
  if (message) {
    printf("Output file(s):\n")
    for (var in names(metGen$settings$outVars)) {
      printf("%15s = %s\n", var, metGen$settings$outVars[[var]]$filename)
    }
  }
}


mgsetInit<- function() {
  mgsetInitSettings()
  mgsetInitMetadata()
  mgsetInitInternal()
}

mgsetInitSettings <- function() {
  metGen$constants <-setConstants()
  metGen$settings <- NULL
  metGen$settings$startDate <- NULL
  metGen$settings$endDate <- NULL
  
  mgsetInDt(24) # Set N hours per timestep
  mgsetOutDt(6) # Set N hours per timestep
  
}

mgsetInitMetadata <- function() {
  assign("metadata", list(), env=metGen)
  metGen$metadata$outvars <- list(
    pr         = list(filename = "", enable = FALSE, units = "mm s-1",    longName = "incoming precipitation"),
    tas        = list(filename = "", enable = FALSE, units = "C",         longName = "air temperature"),
    shortwave  = list(filename = "", enable = FALSE, units = "W m-2",     longName = "incoming shortwave"),
    longwave   = list(filename = "", enable = FALSE, units = "W m-2",     longName = "incoming longwave"),
    pressure   = list(filename = "", enable = FALSE, units = "kPa",       longName = "near surface atmospheric pressure"),
    qair       = list(filename = "", enable = FALSE, units = "kg kg-1",   longName = "specific humidity"),
    vp         = list(filename = "", enable = FALSE, units = "kPa",       longName = "near surface vapor pressure"),
    # relhum     = list(filename = "", enable = FALSE, units = "%",         longName = "relative humidity"),
    relhum     = list(filename = "", enable = FALSE, units = "% / 0.01",         longName = "relative humidity"), ## relhum needs to be fraction. Because fraction does not extist in udunits we call it "% / 0.01"
    density    = list(filename = "", enable = FALSE, units = "kg m-3",    longName = "near-surface atmospheric density"),
    wind       = list(filename = "", enable = FALSE, units = "m s-1",     longName = "near surface wind speed")
  )
  
  metGen$metadata$invars <- list(
    pr         = list(units = "mm s-1",   longName = "incoming precipitation"),
    tasmin     = list(units = "Celsius",  longName = "minimum air temperature"),
    tasmax     = list(units = "Celsius",  longName = "maximum air temperature"),
    shortwave  = list(units = "W m-2",    longName = "shortwave radiation"),
    vp         = list(units = "kPa",      longName = "near surface vapor pressure"),
    relhum     = list(units = "% / 0.01", longName = "relative humidity"), ## relhum needs to be fraction. Because fraction does not extist in udunits we call it "% / 0.01"
    # relhum     = list(units = "fraction", longName = "relative humidity"),
    longwave   = list(units = "W m-2",    longName = "longwave radiation"),
    pressure   = list(units = "kPa",      longName = "near surface atmospheric pressure"),
    wind       = list(units = "m s-1",    longName = "near surface wind speed")
  )
  
  metGen$metadata$elevation <- list(ncName = "elevation")
  # settings$elevation <- list(ncFileName = ncFileNameElevation, ncName = "elevation")
}

mgsetInitInternal <- function() {
  assign("internal", list(), env=metGen)
  
  ## Return the location of the Example NetCDF-files
  metGen$internal$ncFileNameElevation  <- system.file("extdata", "elevation_Mekong.nc4", package = "metGeneratoR")
  metGen$internal$ncFileNamePr         <- system.file("extdata", "pr_Mekong.nc4", package = "metGeneratoR")
  metGen$internal$ncFileNameTasmin     <- system.file("extdata", "tasmin_Mekong.nc4", package = "metGeneratoR")
  metGen$internal$ncFileNameTasmax     <- system.file("extdata", "tasmax_Mekong.nc4", package = "metGeneratoR")
  metGen$internal$ncFileNameWind       <- system.file("extdata", "wind_Mekong.nc4", package = "metGeneratoR")
}

mgsetInit()
