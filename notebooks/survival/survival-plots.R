##################################################
# Survival plotting  #
##################################################

Sys.setenv(LANGUAGE = "en")
install.packages("forestploter", repos="https://cloud.r-project.org")

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
library(forestploter)
library(broom)
library(grid)


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

forest_with_strata <- function(fit,
                               strata_label = "Cancer type",
                               var_labels   = NULL,
                               base_size    = 7,
                               ci_width     = 50) {

  n_strata <- function(fit) {
    mf <- model.frame(fit); scol <- grep("^strata\\(", names(mf), value = TRUE)
    if (length(scol)) return(nlevels(droplevels(as.factor(mf[[scol[1]]]))))
    if (!is.null(fit$strata)) length(fit$strata) else NA_integer_
  }

  mf  <- model.frame(fit)
  td  <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tl  <- attr(fit$terms, "term.labels"); tl <- tl[!grepl("^strata\\(", tl)]

  ## build rows per variable, inserting reference rows for factors
  rows <- lapply(tl, function(v) {
    col <- mf[[v]]
    if (is.factor(col)) {
      levs <- levels(droplevels(col))
      ref  <- levs[1]                                   # reference = first level
      out  <- lapply(levs, function(L) {
        if (L == ref) {
          tibble(Variable = v, Value = L,
                 est = 1, lower = NA, upper = NA, p = NA_real_, is_ref = TRUE)
        } else {
          hit <- td[td$term == paste0(v, L), ]
          tibble(Variable = v, Value = L,
                 est = hit$estimate, lower = hit$conf.low, upper = hit$conf.high,
                 p = hit$p.value, is_ref = FALSE)
        }
      })
      bind_rows(out)
    } else {
      hit <- td[td$term == v, ]                          # continuous: one row, no level
      tibble(Variable = v, Value = "",
             est = hit$estimate, lower = hit$conf.low, upper = hit$conf.high,
             p = hit$p.value, is_ref = FALSE)
    }
  })
  dat <- bind_rows(rows) %>%
    group_by(Variable) %>%
    mutate(Variable = ifelse(row_number() == 1, Variable, "")) %>%
    ungroup()

  if (!is.null(var_labels))
    dat <- dat %>% mutate(Variable = ifelse(Variable == "", "",
                                            recode(Variable, !!!var_labels)))

  n_str <- n_strata(fit)
  dat <- bind_rows(dat, tibble(Variable = strata_label, Value = "",
                               est = NA, lower = NA, upper = NA, p = NA_real_,
                               is_ref = FALSE))

  dat <- dat %>% mutate(
    `HR (95% CI)` = case_when(
      is.na(est) & Variable == strata_label ~ sprintf("stratified (%d strata)", n_str),
      is_ref                                ~ "reference",
      TRUE ~ sprintf("%.2f (%.2f\u2013%.2f)", est, lower, upper)),
    p = ifelse(is.na(p), "",
               paste(ifelse(p < .001, "<0.001", sprintf("%.3f", p)),
                      ifelse(p < .001, "***",
                      ifelse(p < .01,  "**",
                      ifelse(p < .05,  "*", ""))))),
    ` ` = strrep(" ", ci_width))

  tab <- dat %>% select(Variable, Value, ` `, `HR (95% CI)`, p)

  tm <- forest_theme(base_size = base_size,
                     core     = list(bg_params = list(fill = c("grey92","white"))),
                     colhead  = list(fg_params = list(fontface = c(1,1,1,1,3))))  # italic last col header

  forest(tab,
         est = dat$est, lower = dat$lower, upper = dat$upper,
         ci_column = 3, ref_line = 1, x_trans = "log10",
         footnote = sprintf("n = %d, events = %d", fit$n, fit$nevent),
         theme    = tm)
}

save_forestploter <- function(p, file, width = NULL, height = NULL,
                        dpi = 300, scale = 1) {
  # size the device to the plot's natural dimensions unless overridden
  wh <- forestploter::get_wh(plot = p, unit = "in")
  w  <- (if (is.null(width))  wh[1] else width)  * scale
  h  <- (if (is.null(height)) wh[2] else height) * scale

  ext <- tolower(tools::file_ext(file))
  switch(ext,
    svg = svg(file, width = w, height = h),                       # vector; no dpi
    png = png(file, width = w, height = h, units = "in", res = dpi),
    stop("Unsupported extension: ", ext, " (use .svg or .png)"))

  on.exit(dev.off(), add = TRUE)      # guarantees the device closes even on error
  plot(p)                             # forestploter's plot method draws the grob
  invisible(file)
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
    labels = names(survObj$strata) %>% str_replace("^.*=","")
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