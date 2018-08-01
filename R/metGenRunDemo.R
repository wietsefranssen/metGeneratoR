metGenRunDemo <- function() {
  
  mgsetLonlatbox(c(92.25, 110.25, 7.25, 36.25))
  mgsetPeriod(startdate = "1950-1-1", enddate = "1950-1-31")
  mgsetInDt(24) # Set N hours per timestep
  mgsetOutDt(6) # Set N hours per timestep
  
  metGen$constants<-setConstants()
  constants <- metGen$constants
  
  mgsetInVars(list(
    pr         = list(ncname = "pr",      filename = metGen$internal$ncFileNamePr),
    tasmin     = list(ncname = "tasmin",  filename = metGen$internal$ncFileNameTasmin),
    tasmax     = list(ncname = "tasmax",  filename = metGen$internal$ncFileNameTasmax),
    wind       = list(ncname = "wind",      filename = metGen$internal$ncFileNameWind)
  ))
  
  ## Define elevation file
  mgsetElevation(ncname = "elevation", filename = metGen$internal$ncFileNameElevation)
  
  mgsetOutVars(c("tas",
                 "pr",
                 "wind"))
  
  metGenRun()
}
