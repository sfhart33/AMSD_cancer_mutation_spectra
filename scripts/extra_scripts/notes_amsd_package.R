
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
  usethis::use_github_release()
  usethis::use_mit_license("Samuel Hart")
  
# Other things
  usethis::use_testthat()
  usethis::use_test("amsd")
  usethis::use_github_action_check_standard()
  usethis::use_news_md()
  usethis::use_citation()
  usethis::use_version("minor")  # e.g., bump to 0.2.0 before paper release
  usethis::use_git_tag("v0.2.0", message = "Release version 0.2.0")
  usethis::git_push(tag = TRUE)
  
  
  # make a vignette