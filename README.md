# metGeneratoR

## Installation

In Ubuntu/Debian you need to install some prerequisites first (testen on 18.04 LTE):

    sudo apt update
    sudo apt install libssl-dev libcurl4-openssl-dev netcdf-bin libnetcdf-dev udunits-bin libudunits2-dev

Start R and run the following commands:

    library(devtools)
    install_git("https://github.com/wietsefranssen/metGeneratoR")

## Loading

    ## Load the library
    library(metGeneratoR)

## Help

    ## Show help
    ?metGenerator
    
## Usage

Run the following command to test if the package is working properly.

    ## Run the demos from help
    example("metGeneratoR", echo = F)

## LICENSE

MIT
