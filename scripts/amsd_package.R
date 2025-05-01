
setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github")
library(devtools)
library(roxygen2)
devtools::create("amSpecDist")

# Edited DESCRIPTION file
# added functions

# install and test package
  setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/amSpecDist")
  devtools::load_all()


# STILL to Do


  # document functions, then devtools::document()
  devtools::document()
  # make example data available
  # make a vignette
  # put on github