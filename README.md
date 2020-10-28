# metGeneratoR

## Installation

In Ubuntu/Debian you need to install some prerequisites first (testen on 18.04 LTE):

    sudo apt update
    sudo apt install libssl-dev libcurl4-openssl-dev netcdf-bin libnetcdf-dev udunits-bin libudunits2-dev libxml2-dev 
    
Start R and run the following commands:

    devtools::install_git("https://github.com/wietsefranssen/metGeneratoR", ref="waterSIS_2dimlonlat")

## Loading

    library(metGeneratoR)

## Check version

    packageVersion("metGeneratoR")

## Help

    ?metGenerator
    
## Usage

Run the following command to test if the package is working properly.

    example("metGeneratoR", echo = F)

## LICENSE

MIT
