# metGeneratoR

## Installation

Start R and run the following commands:

    library(devtools)
    install_git("https://github.com/wietsefranssen/metGeneratoR")

## Examples

Run the following command to test if the package is working properly.
A file called `example_output_Mekong.nc` should be written in the working directory.

    ## Load the library
    library(metGeneratoR)
    
    ## Run metGen
    metGenRunDemo()

## Usage

TODO: make this better! far from complete!

Below is an example of how you can run the script with settings you define yourself. 

    ## Load the library
    library(metGeneratoR)

    ## Init the settings datatype
    settings <- initSettings(startdate = "1950-01-01",
                             enddate = "1950-1-31",
                             outstep = 6,
                             lonlatbox = c(92.25, 110.25, 7.25, 36.25),
                             outfile = "example_output_Mekong.nc")
  
    ## Some system settings
    settings$system$nCores <- 2
    settings$system$maxMem <- 0.0040 # in Gb
  
    ## Input variables
    ## Comment out the ones you dont want to include
    settings <- setInputVars(settings,list(
      pr         = list(ncFileName = "pr_Mekong.nc4",        ncName = "pr",        vicIndex = 9,   alma = FALSE),
      tasmin     = list(ncFileName = "tasmin_Mekong.nc4",    ncName = "tasmin",    vicIndex = 17),
      tasmax     = list(ncFileName = "tasmax_Mekong.nc4",    ncName = "tasmax",    vicIndex = 16),
      wind       = list(ncFileName = "wind_Mekong.nc4",      ncName = "wind",      vicIndex = 20)
    ))
    settings$elevation <- list(ncFileName = "elevation_Mekong.nc4", ncName = "elevation")
  
    ## Output variables
    ## Comment out the ones you dont want to include
    settings$outputVars <- list(
      pr         = list(VICName = "OUT_PREC",       units = "mm"),
      tas        = list(VICName = "OUT_AIR_TEMP",   units = "C"),
      shortwave  = list(VICName = "OUT_SHORTWAVE",  units = "W m-2"),
      longwave   = list(VICName = "OUT_LONGWAVE",   units = "W m-2"),
      pressure   = list(VICName = "OUT_PRESSURE",   units = "kPa"),
      qair       = list(VICName = "OUT_QAIR",       units = "kg kg-1"),
      vp         = list(VICName = "OUT_VP",         units = "kPa"),
      rel_humid  = list(VICName = "OUT_REL_HUMID",  units = "fraction"),
      density    = list(VICName = "OUT_DENSITY",    units = "kg m-3"),
      wind       = list(VICName = "OUT_WIND",       units = "m s-1")
    )
      
    ## Run the main routine
    mtclim_run(settings)

## LICENSE

TODO
