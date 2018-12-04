#' @export
metGen <- new.env()

#' @export
mgcheckVariables <- function() {
  if (!is.null(metGen$settings$outVars[["vp"]])) {
    if (metGen$metadata$inVars$relhum$enabled && metGen$metadata$inVars$qair$enabled) {
      stop(paste0("Both \"relhum\" and \"qair\" are provided as input. \n",
                  "Select only of the two variables for the calculation of \"vp\""))
    }
    if (!metGen$metadata$inVars$relhum$enabled && !metGen$metadata$inVars$qair$enabled && !metGen$metadata$inVars$vp$enabled ) {
      stop(paste0("\"relhum\" or \"qair\"",
                  " need to be provided as input for the calculation of \"vp\""))
    }
  }
}

#' @export
mgsetLonlatbox <- function(lonlatbox) {
  metGen$settings$lonlatbox <- lonlatbox
  metGen$settings$x <- seq(lonlatbox[1],lonlatbox[2], 0.5)
  metGen$settings$y <- seq(lonlatbox[3],lonlatbox[4], 0.5)
  metGen$settings$nx <- length(metGen$settings$x)
  metGen$settings$ny <- length(metGen$settings$y)
}

#' @export
mgsetInDt <- function(inDt) {
  metGen$derived$inDt <- metGen$settings$inDt <- inDt
  metGen$derived$nInStepDay  <- 24 / metGen$derived$inDt
  
  ## Update mgsetPeriod
  if (!is.null(metGen$settings$startDate) && !is.null(metGen$settings$endDate))
  {
    mgsetPeriod(metGen$settings$startDate, metGen$settings$endDate)
  }
}

#' @export
mgsetOutDt <- function(outDt) {
  metGen$derived$outDt <- metGen$settings$outDt <- outDt
  metGen$derived$nOutStepDay  <- 24 / metGen$derived$outDt
  
  ## Update mgsetPeriod
  if (!is.null(metGen$settings$startDate) && !is.null(metGen$settings$endDate))
  {
    mgsetPeriod(metGen$settings$startDate, metGen$settings$endDate)
  }
}

#' @export
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

#' @export
mgsetElevation <- function(ncname, filename) {
  metGen$settings$elevation$ncname <- ncname
  metGen$settings$elevation$filename <- filename
}

#' @export
mgsetInVars <- function(varlist) {
  metGen$settings$inVars <- NULL
  for (var in names(metGen$metadata$inVars)) {
    metGen$metadata$inVars[[var]]$enabled <- FALSE
  }
  for (var in names(varlist)) {
    ## Check if input var exists
    if (!is.element(var, names(metGen$metadata$inVars))) {
      stop(paste0("The variable ", var, " is not a standard input variable!\n",
                  "Valid variables are:\n\t", 
                  paste(names(metGen$metadata$inVars),collapse=", ")), call. = FALSE)
    } else {
      metGen$metadata$inVars[[var]]$enabled <- TRUE
      metGen$settings$inVars[[var]]$ncname <- varlist[[var]]$ncname
      metGen$settings$inVars[[var]]$filename <- varlist[[var]]$filename
      metGen$settings$inVars[[var]]$longName <- metGen$metadata$inVars[[var]]$longName
    }
  }
  mgcheckInVars()
}

#' @export
mgsetOutVars <- function(varnames) {
  ## Check if inputvars are already defined
  if (length(metGen$settings$inVars) <= 0) {
    cat(paste0("First set the input variables. Based on that the possible output variables can be defined\n"))
    stop()
  } else {
    metGen$settings$outVars <- NULL
    for (var in varnames) {
      ## Check input variables
      if (var == "swdown" && !metGen$metadata$inVars$swdown$enabled) {
        stop(paste0("swdown output can only be generated if swdown is provided as input\nOn the moment, only the following input variables are defined: ", paste(names(metGen$settings$inVars), collapse=", "), "\n"), call. = FALSE)
      }
      if (var == "pr") {
        errmess <- paste0("Precipitation output can only be generated if:\n", 
                          "\tprecipitation is provided as input or\n",
                          "\trainfall and snowfall are provided as input\n",
                          "On the moment, only the following input variables are defined: \n\t", 
                          paste(names(metGen$settings$inVars), collapse=", "), "\n")
        if (metGen$metadata$inVars$pr$enabled) {
          if (metGen$metadata$inVars$rainf$enabled || metGen$metadata$inVars$snowf$enabled) {
            stop(errmess, call. = FALSE)
          }
        } else if (metGen$metadata$inVars$rainf$enabled && metGen$metadata$inVars$snowf$enabled) {
        } else stop(errmess, call. = FALSE)
      }
      
      ## Check if output var exists
      if (!is.element(var, names(metGen$metadata$outVars))) {
        cat(paste0("WARING: The variable ", var, " is not a standard output variable!\n",
                   "        Valid variables are: ", paste(names(metGen$metadata$outVars), collapse=", ") , "\n",
                   "        Ingoring: ", var, "\n"))
      } else {
        metGen$settings$outVars[[var]]$name <- var
        metGen$settings$outVars[[var]]$longName <- metGen$metadata$outVars[[var]]$longName
        metGen$settings$outVars[[var]]$units <- metGen$metadata$outVars[[var]]$output_units
      }
    }
  }
  mgsetOutName(nameString = "output/<VAR>/<VAR>_<SYEAR><SMONTH><SDAY>_<EYEAR><EMONTH><EDAY>.nc", message = F)
  mgcheckOutVars()
  
}

#' @export
mgsetOutName <- function(nameString, message = F) {
  # nameString <- "output/<VAR>/<VAR>_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_landonly_<SYEAR>.nc"
  
  for (var in names(metGen$settings$outVars)) {
    nameStringTmp <- nameString
    ## Substitute Dates
    nameStringTmp <- gsub("<SYEAR>", format(metGen$derived$startDate,format="%Y"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<EYEAR>", format(metGen$derived$endDate,format="%Y"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<SMONTH>", format(metGen$derived$startDate,format="%m"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<EMONTH>", format(metGen$derived$endDate,format="%m"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<SDAY>", format(metGen$derived$startDate,format="%d"), nameStringTmp, ignore.case = T)
    nameStringTmp <- gsub("<EDAY>", format(metGen$derived$endDate,format="%d"), nameStringTmp, ignore.case = T)
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


#'Initialize and set some settings required for the metGen package to run
#'
#' @return nothing dfgdg
#' @examples 
#' mgsetInit()
#' @export
mgsetInit <- function() {
  mgsetInitSettings()
  mgsetInitMetadata()
  mgsetInitInternal()
  mgsetInitConstants()
}

#' @export
mgsetInitSettings <- function() {
  # metGen$constants <- setConstants()
  metGen$settings <- NULL
  metGen$settings$startDate <- NULL
  metGen$settings$endDate <- NULL
  
  mgsetInDt(24) # Set N hours per timestep
  mgsetOutDt(6) # Set N hours per timestep
}

#' @export
mgsetInitMetadata <- function() {
  assign("metadata", list(), env=metGen)
  
  ## Needed input unit
  metGen$metadata$inVars <- list(
    pr         = list(input_units = "", internal_units = "mm s-1"),   # incoming precipitation 
    rainf      = list(input_units = "", internal_units = "mm s-1"),   # incoming rainfall 
    snowf      = list(input_units = "", internal_units = "mm s-1"),   # incoming snowfall 
    tas        = list(input_units = "", internal_units  = "Celsius"),  # average air temperature 
    tasmin     = list(input_units = "", internal_units  = "Celsius"),  # minimum air temperature 
    tasmax     = list(input_units = "", internal_units  = "Celsius"),  # maximum air temperature 
    swdown     = list(input_units = "", internal_units  = "W m-2"),    # shortwave radiation 
    vp         = list(input_units = "", internal_units  = "kPa"),      # near surface vapor pressure 
    relhum     = list(input_units = "", internal_units  = "% / 0.01"), # relative humidity  ## relhum needs to be fraction. Because fraction does not extist in udunits  we call it "% / 0.01"
    # relhum     = list(input_units = "", internal_units  = "fraction"), # relative humidity 
    qair       = list(input_units = "", internal_units  = "kg/kg"),   # near surface specific humidity 
    lwdown     = list(input_units = "", internal_units  = "W m-2"),   # longwave radiation 
    psurf   = list(input_units = "", internal_units  = "kPa"),     # near surface atmospheric pressure 
    wind       = list(input_units = "", internal_units  = "m s-1")    # near surface wind speed
  )
  
  ## Output metadata (output_units)
  metGen$metadata$outVars <- list(
    pr         = list(filename = "", enable = FALSE, internal_units = "mm s-1",    output_units = "mm",        longName = "incoming precipitation"),
    tas        = list(filename = "", enable = FALSE, internal_units = "C",         output_units = "C",         longName = "air temperature"),
    tasmin     = list(filename = "", enable = FALSE, internal_units = "C",         output_units = "C",         longName = "minimum air temperature"),
    tasmax     = list(filename = "", enable = FALSE, internal_units = "C",         output_units = "C",         longName = "maximum air temperature"),
    swdown     = list(filename = "", enable = FALSE, internal_units = "W m-2",     output_units = "W m-2",     longName = "incoming shortwave"),
    lwdown     = list(filename = "", enable = FALSE, internal_units = "W m-2",     output_units = "W m-2",     longName = "incoming longwave"),
    psurf      = list(filename = "", enable = FALSE, internal_units = "kPa",       output_units = "kPa",       longName = "near surface atmospheric pressure"),
    qair       = list(filename = "", enable = FALSE, internal_units = "kg kg-1",   output_units = "kg kg-1",   longName = "specific humidity"),
    vp         = list(filename = "", enable = FALSE, internal_units = "kPa",       output_units = "kPa",       longName = "near surface vapor pressure"),
    # relhum     = list(filename = "", enable = FALSE, internal_units = "%",         output_units = "%",         longName = "relative humidity"),
    relhum     = list(filename = "", enable = FALSE, internal_units = "% / 0.01",  output_units = "% / 0.01",  longName = "relative humidity"), ## relhum needs to be fraction. Because fraction does not extist in udunits we call it "% / 0.01"
    density    = list(filename = "", enable = FALSE, internal_units = "kg m-3",    output_units = "kg m-3",    longName = "near-surface atmospheric density"),
    wind       = list(filename = "", enable = FALSE, internal_units = "m s-1",     output_units = "m s-1",     longName = "near surface wind speed")
  )
  
  metGen$metadata$elevation <- list(ncName = "elevation")
  # settings$elevation <- list(ncFileName = ncFileNameElevation, ncName = "elevation")
}

#' @export
mgsetInitInternal <- function() {
  assign("internal", list(), env=metGen)
  
  ## Return the location of the Example NetCDF-files
  metGen$internal$ncFileNamePr         <- system.file("extdata", "pr_19500101_19500131.nc4",      package = "metGeneratoR")
  metGen$internal$ncFileNameTasmin     <- system.file("extdata", "tasmin_19500101_19500131.nc4",  package = "metGeneratoR")
  metGen$internal$ncFileNameTasmax     <- system.file("extdata", "tasmax_19500101_19500131.nc4",  package = "metGeneratoR")
  metGen$internal$ncFileNamePs         <- system.file("extdata", "ps_19500101_19500131.nc4",      package = "metGeneratoR")
  metGen$internal$ncFileNameswdown     <- system.file("extdata", "rsds_19500101_19500131.nc4",    package = "metGeneratoR")
  metGen$internal$ncFileNamelwdown     <- system.file("extdata", "rlds_19500101_19500131.nc4",    package = "metGeneratoR")
  metGen$internal$ncFileNameRelhum     <- system.file("extdata", "hurs_19500101_19500131.nc4",    package = "metGeneratoR")
  metGen$internal$ncFileNameWind       <- system.file("extdata", "sfcWind_19500101_19500131.nc4", package = "metGeneratoR")
}

#' @export
mgsetInitConstants <- function() {
  metGen$constants <- setConstants()
}

