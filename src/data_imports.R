##################################################
# Supplementary Table imports in R.              #
# See data_imports.py for the python equivalents.#
##################################################

Sys.setenv(LANGUAGE = "en")

library(tidyverse)
library(readxl)

patients <- function(path){
    ## path: path to data/Supplementary Tables.xlsx
    ## or data/Supplementary Tables 12_1_24.xlsx
    tbl <- read_excel(path, sheet="1. Patients",
                     col_types=c(rep('text',4),'numeric',rep('text',5),'numeric')) %>% 
        suppressWarnings()
    return(tbl)
}