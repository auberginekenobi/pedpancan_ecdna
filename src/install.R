# Install all dependencies necessary to run this project.

install.packages("renv")
library(renv)
renv::init(bare = TRUE)
renv::activate()

# Required packages
install.packages("pacman")
pacman::p_load(
  
  # project and file management
  #############################
  here,     # file paths relative to R project root folder
  rio,      # import/export of many types of data

  # package install and management
  ################################
  pacman,   # package install/load
  renv,     # managing versions of packages when working in collaborative groups
  remotes,  # install from github
  yaml,
  
  # General data management
  #########################
  tidyverse)   # includes many packages for tidy data wrangling and presentation
pacman::p_load(  
  naniar,      # assessing missing data
  labelled,    # 
  writexl,
  
  # statistics  
  ############
  janitor,      # tables and data cleaning
  gt,
  gtsummary,    # making descriptive and statistical tables
  survival,     # KM and Cox models

  # plots - general
  #################
  RColorBrewer,     # color scales
  survminer,        # ggplot2 survival curves
  ggsurvfit,
  partykit,         # visualizing tree-structured regression and classification models.
  rmarkdown,
  cowplot,
  extrafont
)

renv::snapshot()
