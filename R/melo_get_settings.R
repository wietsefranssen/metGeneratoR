melo_get_settings <- function() {
  if (!interactive()) {
    args = commandArgs(trailingOnly=TRUE)
    
    # Test arguments
    if (length(args)<6) {
      stop("Too less arguments. \n\nRequired:\n  --> varType, projection, nhourly, inFile, varname, outfile", call.=FALSE)
    } else if (args[1]=="tair") {
      if (length(args)<8) {
        stop("Too less arguments. \n\nRequired:\n  --> varType, projection, nhourly, inFileTmin, varname, inFileTmax, varname, outfile\n\n ex: ./melo_run.R tair xy 3 ~/tn_19900101.nc tn ~/tx_19900101.nc tx ~/out.nc", call.=FALSE)
      }
    } else if (length(args)==2) {
      # default output file
      print(args[2])
    }
    
    varType <- args[1]
    fileType <- args[2]
    nhourly <- as.numeric(args[3])
    inFile <- args[4]
    varname <- args[5]
    if (varType == "tair") {
      inFile2 <- args[6]
      varname2 <- args[7]
      outFile <- args[8]
    } else {
      outFile <- args[6]
    }
  } else {
    # varType <- "tair"
    varType <- "potrad"
    fileType <- "xy"
    nhourly <- 1
    inFile <-"~/tn_19900101.nc"
    varname <- "tn"
    
    if (varType == "tair") {
      inFile2 <- "~/tx_19900101.nc"
      varname2 <- "tx"
    }
    outFile <- "~/out_interactive.nc"
  }
  
  ## Check
  if (varType != "potrad" && varType != "precip" && varType != "tair" && varType != "shortwave") stop("varType not supported!")
  
  if (varType == "potrad" || varType == "tair" ) {
    outVar <- varType
  } else {
    outVar <- varname
  }
  
  ## Print setting
  cat(paste("\n######################################################################################"))
  cat(paste("\n## varType:   ", varType, "    \tfileType:   ", fileType, "   nhourly:    ", nhourly))
  cat(paste("\n## inFile:   ", inFile,"    varname:   ", varname))
  if (varType == "tair") {
    cat(paste("\n## inFile2:   ", inFile2,"    varname2:   ", varname2))
  }
  cat(paste("\n## outFile:   ", outFile))
  cat(paste("      outVar:   ", outVar))
  cat(paste("\n######################################################################################\n"))
  
  ## Create function output
  settings <- NULL
  settings$fileType <- fileType
  settings$outFile <- outFile
  settings$varType <- varType
  settings$nhourly <- nhourly
  settings$inFile <- inFile
  settings$varname <- varname
  settings$outVar <- outVar
  if (varType == "tair") settings$inFile2 <- inFile2
  if (varType == "tair") settings$varname2 <- varname2
  
  return(settings)
}
