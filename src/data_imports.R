##################################################
# Supplementary Table imports in R.              #
# See data_imports.py for the python equivalents.#
##################################################

Sys.setenv(LANGUAGE = "en")

library(tidyverse)
library(readxl)

patients <- function(path){
    ## path: path to data/Supplementary Tables.xlsx
    tbl <- read_excel(path, sheet="1. Patients")
    return(tbl)
}
biosamples <- function(path){
    ## path: path to data/Supplementary Tables.xlsx
    tbl <- read_excel(path, sheet="2. Biosamples")
    return(tbl)
}
genes <- function(path){
    ## path: path to data/Supplementary Tables.xlsx
    tbl <- read_excel(path, sheet="5. Gene amplifications", col_types="text")
    return(tbl)
}