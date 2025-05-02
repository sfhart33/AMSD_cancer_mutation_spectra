
setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github")
library(devtools)
library(roxygen2)
library(usethis)


devtools::create("mutspecdist")

# Edited DESCRIPTION file
# added functions

devtools::document()

# install and test package
  setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/mutspecdist")
  devtools::load_all()

# make example data available
  mouse_sample_table <- sample_table
  usethis::use_data(mouse_carcinogen_spectra)
  usethis::use_data(mouse_sample_table)
  
# put on github
  usethis::use_git()
  usethis::use_github()

  
  
  
  
  

  
  
  
  # make a vignette