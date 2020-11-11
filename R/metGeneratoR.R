#' metGeneratoR package
#' 
#' The metGeneratoR package disaggregates daily NetCDF data to subdaily data.
#' This can be 1hourly, 2hourly, 3hourly, 4hourly, 6hourly or 8hourly.
#' If needed it converts data to the output unit needed for VIC.
#' 
#' @docType package
#' @keywords internal
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib metGeneratoR
#' @name metGeneratoR
#' @examples
#' ########################
#' ## preprocessing
#' ## TODO
#' ########################
#' 
#' ## (re-)Initialize the metGen environment with basic settings
#' ##
#' ## Please always rerun this function before making adjustments in the metGeneratoR setup
#' mgsetInit()
#' 
#' ## Set the period to run
#' mgsetPeriod(startdate = "1950-1-1", enddate = "1950-1-05")
#' 
#' ## Set the output timesteps in hours
#' mgsetOutDt(6)  # 6 hourly output
#' 
#' ## Set the input variables
#' mgsetInVars(list(
#'   pr         = list(ncname = "pr",      filename = metGen$internal$ncFileNamePr)
#' ))
#' 
#' ## Set the output variables
#' mgsetOutVars(c("radfrac", "tminhour", "tmaxhour"))
#' 
#' ## Run the metGeneratoR based on the settings above
#' ## The output files will be written to (a subfolder of) the current working directory of R
#' metGenRun()
#' 
#' ########################
#' ## do the job (END)  
#' ## 
#' ## The example below 
#' ########################
#' 
#' ## (re-)Initialize the metGen environment with basic settings
#' ##
#' ## Please always rerun this function before making adjustments in the metGeneratoR setup
#' mgsetInit()
#' 
#' ## Set the period to run
#' mgsetPeriod(startdate = "1950-1-1", enddate = "1950-1-05")
#' 
#' ## Set the output timesteps in hours
#' mgsetOutDt(6)  # 6 hourly output
#' 
#' ## Set the input variables
#' mgsetInVars(list(
#'   pr         = list(ncname = "pr",      filename = metGen$internal$ncFileNamePr),
#'   tasmin     = list(ncname = "tasmin",  filename = metGen$internal$ncFileNameTasmin),
#'   tasmax     = list(ncname = "tasmax",  filename = metGen$internal$ncFileNameTasmax),
#'   wind       = list(ncname = "sfcWind", filename = metGen$internal$ncFileNameWind),
#'   swdown     = list(ncname = "rsds",    filename = metGen$internal$ncFileNameswdown),
#'   lwdown     = list(ncname = "rlds",    filename = metGen$internal$ncFileNamelwdown),
#'   relhum     = list(ncname = "hurs",    filename = metGen$internal$ncFileNameRelhum),
#'   psurf      = list(ncname = "ps",      filename = metGen$internal$ncFileNamePs),
#'   radfrac    = list(ncname = "radfrac",  filename = "./radfrac_19500101_19500105.nc"),
#'   tminhour   = list(ncname = "tminhour", filename = "./tminhour_19500101_19500105.nc"),
#'   tmaxhour   = list(ncname = "tmaxhour", filename = "./tmaxhour_19500101_19500105.nc")
#' ))
#' 
#' ## Set the output variables
#' mgsetOutVars(c("tas", "pr", "wind", "vp", "psurf", "swdown", "lwdown"))
#' 
#' ## Run the metGeneratoR based on the settings above
#' ## The output files will be written to (a subfolder of) the current working directory of R
#' metGenRun()
mgsetStart()
