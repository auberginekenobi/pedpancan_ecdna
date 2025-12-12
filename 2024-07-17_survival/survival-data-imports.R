##################################################
# Survival data loading and preprocessing #
##################################################

Sys.setenv(LANGUAGE = "en")

library(tidyverse)
library(readxl)
library(dplyr)
library(stringr)
library(naniar) #for replace with Nas function

load_suppl_tbl_1 <- function(path){
    ## path: path to data/Supplementary Tables.xlsx
    ## or data/Supplementary Tables 12_1_24.xlsx
    tbl <- read_excel(path, sheet="1. Patients") %>%
        mutate(OS_months = suppressWarnings(as.numeric(OS_months)),
               age_at_diagnosis = suppressWarnings(as.numeric(age_at_diagnosis)))
    return(tbl)
}
load_chapman_2023 <- function(path,include_archer=TRUE){
    ## path: path to Chapman et al. 2023 supplementary tables.
    ## include_archer: whether to include the "Archer" subcohort. These should be included (they're ICGC samples), but were not until 2025-11-06.
    
    tbl <- read_excel(path, sheet="1 WGS Patient Cohort")
    
    # unify column names
    tbl = tbl %>% rename(patient_id=Patient_ID, sex=Sex, cancer_subclass=Subgroup, OS_status=Vital_status, cohort=Source, age_at_diagnosis=Age_at_diagnosis)
    
    # get amplification information
    amps <- read_excel(path, sheet="3 Amplicon Classifications") %>%
        mutate(
            amplicon_decomposition_class = case_when(
                str_detect(Notes, regex("blacklist", ignore_case = TRUE)) ~ "No amp/Invalid",
                TRUE ~ amplicon_decomposition_class
            )
        )
    annotate_amplicons <- function(patients, amplicons){
        df <- amplicons %>% 
            rename(patient_id=Patient_ID) %>%
            group_by(patient_id) %>%
            summarize(
                amplicon_class = case_when(
                    any(ecDNA > 0) ~ "ecDNA",
                    any(amplicon_decomposition_class %in% c('Cyclic','Complex non-cyclic','Linear amplification')) ~ "intrachromosomal",
                    TRUE ~ "no amplification"
                ),
                .groups = "drop"
            )
        patients = patients %>%
            left_join(df, by='patient_id') %>%
            mutate(amplicon_class = replace_na(amplicon_class, "no amplification"))
        return(patients)
    }
    tbl = annotate_amplicons(tbl,amps)

    # subset ICGC samples
    if (include_archer) {
        format_archer_ids <- function(id_vector){
            # format Archer IDs to match ICGC
            prefixes <- c("ICGC_", "MBRep", "MDT-AP")
            pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
            case_when( 
                str_detect(id_vector, pattern) ~ id_vector,  # already has prefix → leave unchanged
                TRUE ~ paste0(
                    "ICGC_",
                    str_replace(id_vector, "^([A-Za-z]+)0*([0-9]+)$", "\\1\\2")  # remove leading zeros in number
                )
            )
        }
        tbl = tbl %>% 
            filter(cohort %in% c("ICGC","Archer")) %>%
            mutate(patient_id = format_archer_ids(patient_id))
    } else {
        tbl = tbl %>% filter(cohort == 'ICGC')
    }

    # Standardize units
    tbl = tbl %>%
        mutate(OS_months = Survival_time_years*12) %>%
        mutate(cancer_type = 'MBL') %>%
        mutate(age_at_diagnosis = round(age_at_diagnosis*365.25)) %>%
        mutate(sex=recode(sex, m='Male', f='Female')) %>%
        mutate(OS_status=recode(OS_status, alive='Alive', deceased='Deceased')) %>%
        mutate(cohort = 'ICGC') %>%
        select(c('patient_id','sex','age_at_diagnosis','cohort','cancer_type','cancer_subclass','amplicon_class','OS_status','OS_months'))
    return(tbl)
}

preprocess_survival_data <- function(combinedsurv){
  # Drop NAs
  combinedsurv <- combinedsurv %>%
    filter(complete.cases(amplicon_class,OS_status,OS_months)) %>%
    mutate(OS_months = as.numeric(OS_months)) %>%
  # Censor at 5 years = 60 months
    mutate(OS_months_5y = if_else(OS_months < 60, OS_months, 60)) %>%
    mutate(OS_status_5y = if_else(OS_months <= 60, OS_status, "Alive")) %>%
    mutate(OS_status_5y = if_else(OS_status_5y == "Alive", 0, 1)) %>%
  # get ecDNA status
    mutate(ecDNA_status = if_else(amplicon_class == "ecDNA", "ecDNA+", "ecDNA-")) %>%
    mutate(amplicon_class = if_else(amplicon_class == "intrachromosomal", "chromosomal", amplicon_class)) %>%
    mutate(amplified = if_else(amplicon_class %in% c("ecDNA","chromosomal"), TRUE, FALSE)) %>%
  # zscore age
    mutate(age_at_diagnosis = as.numeric(scale(age_at_diagnosis))) %>%
  # convert to factors
    mutate(ecDNA_status = factor(ecDNA_status)) %>%
    mutate(amplicon_class = factor(amplicon_class)) %>%
    mutate(cancer_type = factor(cancer_type)) %>%
    mutate(amplified = factor(amplified))
  combinedsurv$amplified = relevel(combinedsurv$amplified,ref=TRUE)
    
  return(combinedsurv)
}

load_survival_data <- function(path, mb_path=NULL, include_archer=TRUE){
    # with 1 arg, this function retains the old functionality as of 12/24 and expects the file data/Supplementary Tables 12_1_24.xlsx.
    # with 2 args, this function dynamically builds the same table from data/Supplementary Tables.xlsx and the Chapman et al., 2023 supplementary tables.
    if (is.null(mb_path)) {
        combinedsurv = load_suppl_tbl_1(path)
    } else {
        combinedsurv = bind_rows(
            load_suppl_tbl_1(path),
            load_chapman_2023(mb_path,include_archer)
        )
    }
    return(preprocess_survival_data(combinedsurv))
}

test_load_chapman_2023 <- function(){
    df1 = load_chapman_2023('../data/external/Chapman2023/41588_2023_1551_MOESM4_ESM.xlsx')
    df2 <- load_survival_data("../data/Supplementary Tables 12_1_24.xlsx")
    return(df1 %>%
        inner_join(df2 %>% select(patient_id, class_truth = amplicon_class), by = "patient_id") %>%
        filter(amplicon_class != class_truth)
    )
}