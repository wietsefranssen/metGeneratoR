metGen <- new.env()

mgsetLonlatbox <- function(lonlatbox) {
  metGen$settings$lonlatbox <- lonlatbox
}

mgsetNHourPerStep <- function(nHourPerStep) {
  metGen$settings$nHourPerStep <- nHourPerStep
  metGen$settings$nOutStepDay  <- 24 / nHourPerStep
}

mgsetPeriod <- function(startdate, enddate) {
  metGen$settings$startDate <- startdate
  metGen$settings$endDate  <- enddate
}

mgsetNCores <- function(nCores) {
  metGen$settings$nCores <- nCores
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
      metGen$settings$inVar[[var]]$units <- metGen$metadata$invars[[var]]$units
    }
  }
}

mgsetOutVars <- function(varnames) {
  metGen$settings$outVars <- NULL
  for (var in varnames) {
    ## Check if output var exists
    if (!is.element(var, names(metGen$metadata$outvars))) {
      cat(paste0("WARING: The variable ", var, " is not a standard output variable!\n",
                 "        Valid variables are: ", paste(names(metGen$metadata$outvars),collapse=", ") , "\n",
                 "        Ingoring: ", var, "\n"))
    } else {
      # # if (settings$outperyear) {
      metGen$settings$outVars[[var]]$name <- var
      metGen$settings$outVars[[var]]$filename <- paste0(var, ".nc")
      # metGen$settings$outVars[[var]]$filename <- paste0(var, "_", as.character(settings$ncOut$year), ".nc")
      metGen$settings$outVars[[var]]$longName <- metGen$metadata$outvars[[var]]$longName
      metGen$settings$outVars[[var]]$units <- metGen$metadata$outvars[[var]]$units
      # # } else {
      # # }
    }
  }
}

mgsetInit<- function() {
  mgsetInitSettings()
  mgsetInitMetadata()
  mgsetInitInternal()
}

mgsetInitSettings <- function() {
  metGen$settings <- NULL
  # assign("settings", list(), env=metGen)
  # unlockBinding("settings", metGen)
  mgsetLonlatbox(c(-179.75, 179.75, -89.75, 89.75))
  mgsetNHourPerStep(24) # Set N hours per timestep
  mgsetNCores(1) # Set N hours per timestep
  
}

mgsetInitMetadata <- function() {
  assign("metadata", list(), env=metGen)
  metGen$metadata$outvars <- list(
    pr         = list(filename = "", enable = FALSE, units = "mm",        longName = "incoming precipitation"),
    tas        = list(filename = "", enable = FALSE, units = "C",         longName = "air temperature"),
    shortwave  = list(filename = "", enable = FALSE, units = "W m-2",     longName = "incoming shortwave"),
    longwave   = list(filename = "", enable = FALSE, units = "W m-2",     longName = "incoming longwave"),
    pressure   = list(filename = "", enable = FALSE, units = "kPa",       longName = "near surface atmospheric pressure"),
    qair       = list(filename = "", enable = FALSE, units = "kg kg-1",   longName = "specific humidity"),
    vp         = list(filename = "", enable = FALSE, units = "kPa",       longName = "near surface vapor pressure"),
    rel_humid  = list(filename = "", enable = FALSE, units = "fraction",  longName = "relative humidity"),
    density    = list(filename = "", enable = FALSE, units = "kg m-3",    longName = "near-surface atmospheric density"),
    wind       = list(filename = "", enable = FALSE, units = "m s-1",     longName = "near surface wind speed")
  )
  
  metGen$metadata$invars <- list(
    pr         = list(units = "mm",    longName = "incoming precipitation"),
    tasmin     = list(units = "C",     longName = "minimum air temperature"),
    tasmax     = list(units = "C",     longName = "maximum air temperature"),
    wind       = list(units = "m s-1", longName = "near surface wind speed")
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

# # Gets the value of the variable
# mggetArea <- function() {
#   return(get("area", metGen))
# }
