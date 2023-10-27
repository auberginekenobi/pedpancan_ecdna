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
  ##openxlsx, # import/export of multi-sheet Excel workbooks 

  # package install and management
  ################################
  pacman,   # package install/load
  renv,     # managing versions of packages when working in collaborative groups
  remotes,  # install from github
  yaml,
  
  # General data management
  #########################
  tidyverse,   # includes many packages for tidy data wrangling and presentation
  #dplyr,      # data management
  #tidyr,      # data management
  #ggplot2,    # data visualization
  #stringr,    # work with strings and characters
  #forcats,    # work with factors 
  #lubridate,  # work with dates
  #purrr       # iteration and working with lists
  #readxl      # import/export of Excel workbooks
  naniar,      # assessing missing data
  labelled,    # 
  writexl,
  
  # statistics  
  ############
  janitor,      # tables and data cleaning
  gt,
  gtsummary,    # making descriptive and statistical tables
  ##rstatix,      # quickly run statistical tests and summaries
  ##broom,        # tidy up results from regressions
  ##lmtest,       # likelihood-ratio tests
  ##easystats,
  survival,     # KM and Cox models

  # plots - general
  #################
  #ggplot2,         # included in tidyverse
  RColorBrewer,     # color scales
  survminer,        # ggplot2 survival curves
  ggsurvfit,
  partykit,         # visualizing tree-structured regression and classification models.
  rmarkdown,
)

renv::snapshot()
