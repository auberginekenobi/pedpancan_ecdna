##################################################
# Survival plotting  #
##################################################

Sys.setenv(LANGUAGE = "en")

library(tidyverse)
library(readxl)
library(dplyr)
library(stringr)
library(naniar) #for replace with Nas function
library(survival)
library(survminer)
library(RColorBrewer)
library(janitor)
library(gt)
library(gtsummary)
library(ggsurvfit)
library(extrafont)
library(svglite)

# TODO refactor into separate file
cox_plot <- function(coxobj,data,outfile=NULL,width=3,height=6){
  ## perform a Cox regression and generate the plot
  #coxph(Surv(OS_months, OS_status) ~ ecDNA_status + strata(cancer_type), data = data)
  zph <-cox.zph(coxobj) 
  print(zph)
  #ggcoxzph(zph)
  #m4
  #creating forest plots
  plt <- ggforest(coxobj,data=model.frame(coxobj)) +
        theme_classic(base_size=7, base_family="Arial") +
        theme(axis.text = element_text(size=7,colour="black"),
              plot.title = element_text(size=7))

  if(!is.null(outfile)){
    pngName = paste(outfile, ".png", sep="")
    svgName = paste(outfile, ".svg", sep = "")
    ggsave(path="out", device="png", filename=pngName, width=width, height=height, units='in')
    ggsave(path="out", device="svg", filename=svgName, width=width, height=height, units='in')
  }
  return(plt)
}

km_plot <- function(survObj,outfile=NULL){
  ## perform a KM analysis and generate the plot
  if (length(survObj$n) == 2){
    colors = c('blue', 'red')
    labels = c('ecDNA-', 'ecDNA+')
  } else if (length(survObj$n) == 3){
    colors = c('magenta','red','dodgerblue')
    labels = c('chromosomal','ecDNA','no amplification')
  } else if (length(survObj$n) == 4){
    colors = c('red4','indianred1','orchid4','orchid1')
    colors = c('orchid4','orchid1','red4','indianred1')
    labels = names(survObj$strata) %>% str_replace("^group=","")
  } else {
    stop('colors, labels not defined for this case.')
  }
  plt <- survObj %>% 
   ggsurvfit(linewidth=0.5) +
   labs(x = 'Follow-up time (Months)',
        y = 'Overall Survival') +
   scale_color_manual(values = colors,
                      labels = labels) +
   scale_fill_manual(values = colors,
                     labels = labels) +
   scale_y_continuous(limits=c(0, 1))+
   add_censor_mark(size = .5, alpha = 1) +
   add_risktable(risktable_stats = "n.risk", size=2,
                 theme = theme_risktable_default(axis.text.y.size = 7,
                                    plot.title.size = 7)) +
   add_risktable_strata_symbol(size=4) + 
   theme_classic(base_size=7, base_family="Arial",) +
   theme(axis.text = element_text(size=7,colour="black"),
         legend.position = "bottom",
   )
  if (length(survObj$n) <=3){
    plt <- plt + add_confidence_interval()
  }
  
  if(!is.null(outfile)){
    pngName = paste(outfile, ".png", sep="")
    svgName = paste(outfile, ".svg", sep = "")
    ggsave(path="out", device="png", filename=pngName, width=3, height=3.5, units='in')
    ggsave(path="out", device="svg", filename=svgName, width=3, height=3.5, units='in')
  }
  return(plt)
}