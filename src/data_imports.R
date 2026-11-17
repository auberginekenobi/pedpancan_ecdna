##################################################
# Supplementary Table imports in R.              #
# See data_imports.py for the python equivalents.#
##################################################

Sys.setenv(LANGUAGE = "en")

library(tidyverse)
library(readxl)

bt_types = c('ETMR','GNT','LGG','MNG','HGG','EPN','MBL','CPG','EMBT','ATRT','CHDM','CPT','PINT','CNSTUD','BTNEC')


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

validation_biosamples <- function(path){
    ## path: path to data/Supplementary Tables.xlsx
    tbl <- read_excel(path, sheet="7. Additional tumors")
    return(tbl)
}

amplicon_class_priority <- c('ecDNA', 'chromosomal', 'no amplification') # highest priority first

aggregate_amplicon_class <- function(classes, priority = amplicon_class_priority){
    ## Collapse amplicon classes across biosamples: ecDNA if any biosample has ecDNA; chromosomal if none has ecDNA but
    ## one has chromosomal amplification; otherwise no amplification. Missing values are ignored (NA if all are missing).
    ## classes: character vector of amplicon classes.
    ## priority: classes in descending priority.
    classes <- unique(classes[!is.na(classes)])
    unknown <- setdiff(classes, priority)
    if (length(unknown) > 0) stop('Unknown amplicon classes: ', paste(unknown, collapse = ', '))
    hit <- priority[priority %in% classes]
    if (length(hit) == 0) NA_character_ else hit[1]
}

validation_patients <- function(path){
    ## path: path to data/Supplementary Tables.xlsx
    ## Validation cohorts (ETMR and MYCN-amplified EPN; S.T. 7), one row per patient.
    ## amplicon_class and C19MC_amplicon_class are aggregated over the patient's biosamples with aggregate_amplicon_class;
    ## TTYH1-C19MC_fusion_confirmed is 1 if confirmed in any biosample.
    tbl <- validation_biosamples(path)
    # diagnosis must be the same across a patient's biosamples
    conflicts <- tbl %>%
        group_by(patient_id) %>%
        summarise(across(c(cancer_type, cancer_subclass), ~ n_distinct(.x, na.rm = TRUE)), .groups = 'drop') %>%
        filter(cancer_type > 1 | cancer_subclass > 1)
    if (nrow(conflicts) > 0) stop('Conflicting diagnoses among biosamples from patients: ', paste(conflicts$patient_id, collapse = ', '))
    tbl %>%
        group_by(patient_id) %>%
        summarise(
            cancer_type = first(na.omit(cancer_type)),
            cancer_subclass = first(na.omit(cancer_subclass)),
            amplicon_class = aggregate_amplicon_class(amplicon_class),
            C19MC_amplicon_class = aggregate_amplicon_class(C19MC_amplicon_class),
            `TTYH1-C19MC_fusion_confirmed` = if (all(is.na(`TTYH1-C19MC_fusion_confirmed`))) NA_real_
                                             else max(`TTYH1-C19MC_fusion_confirmed`, na.rm = TRUE),
            n_biosamples = n(),
            .groups = 'drop')
}
