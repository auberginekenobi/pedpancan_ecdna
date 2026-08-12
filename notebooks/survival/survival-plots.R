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


cox_plot <- function(coxobj,data,width=3,height=6){
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
                     colhead  = list(fg_params = list(fontface = c(2,2,2,2,4))))  # bold headers; bold-italic last (p) col

  forest(tab,
         est = dat$est, lower = dat$lower, upper = dat$upper,
         ci_column = 3, ref_line = 1, x_trans = "log10",
         footnote = sprintf("n = %d, events = %d", fit$n, fit$nevent),
         theme    = tm)
}

forest_coxph <- function(fit, var_labels = NULL,
                         n = fit$n, events = fit$nevent,
                         base_size = 7, ci_width = 50) {
  ## Forest plot for a coxph fit with no strata: a reference row for each factor and one row
  ## per continuous term. `n`/`events` feed the footnote and default to the fit's own counts;
  ## override `n` for a person-time-split (counting-process) model, whose fit$n counts
  ## intervals rather than patients. `var_labels` renames the term rows (named vector).
  mf <- model.frame(fit)
  td <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tl <- attr(fit$terms, "term.labels"); tl <- tl[!grepl("^strata\\(", tl)]

  rows <- lapply(tl, function(v) {
    col <- mf[[v]]
    if (is.factor(col)) {
      levs <- levels(droplevels(col)); ref <- levs[1]
      bind_rows(lapply(levs, function(L) {
        if (L == ref) {
          tibble(Variable = v, Value = L, est = 1, lower = NA, upper = NA, p = NA_real_, is_ref = TRUE)
        } else {
          hit <- td[td$term == paste0(v, L), ]
          tibble(Variable = v, Value = L, est = hit$estimate, lower = hit$conf.low,
                 upper = hit$conf.high, p = hit$p.value, is_ref = FALSE)
        }
      }))
    } else {
      hit <- td[td$term == v, ]                          # continuous: one row, no level
      tibble(Variable = v, Value = "", est = hit$estimate, lower = hit$conf.low,
             upper = hit$conf.high, p = hit$p.value, is_ref = FALSE)
    }
  })
  dat <- bind_rows(rows) %>%
    group_by(Variable) %>%
    mutate(Variable = ifelse(row_number() == 1, Variable, "")) %>%
    ungroup()
  if (!is.null(var_labels))
    dat <- dat %>% mutate(Variable = ifelse(Variable == "", "",
                                            recode(Variable, !!!var_labels)))

  dat <- dat %>% mutate(
    `HR (95% CI)` = case_when(
      is_ref ~ "reference",
      TRUE   ~ sprintf("%.2f (%.2f–%.2f)", est, lower, upper)),
    p = ifelse(is.na(p), "",
               paste(ifelse(p < .001, "<0.001", sprintf("%.3f", p)),
                      ifelse(p < .001, "***",
                      ifelse(p < .01,  "**",
                      ifelse(p < .05,  "*", ""))))),
    ` ` = strrep(" ", ci_width))

  tab <- dat %>% select(Variable, Value, ` `, `HR (95% CI)`, p)

  tm <- forest_theme(base_size = base_size,
                     core    = list(bg_params = list(fill = c("grey92","white"))),
                     colhead = list(fg_params = list(fontface = c(2,2,2,2,4))))  # bold headers; bold-italic p col

  forest(tab,
         est = dat$est, lower = dat$lower, upper = dat$upper,
         ci_column = 3, ref_line = 1, x_trans = "log10",
         footnote = sprintf("n = %d, events = %d", n, events),
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

caterpillar_table <- function(coxme_fit, mle_fit, data, focus,
                              group = "molecular_subtype", event = "OS_status_5y") {
  ## Assemble the per-subtype MLE-vs-MAP table consumed by caterpillar_forest(), for any
  ## 0/1 exposure `focus` (e.g. "ecDNA_num" or "amp_num").
  ##   coxme_fit  a coxme with a random slope (0 + focus | group): supplies the pooled
  ##              effect (fixef), between-group variance tau^2 (VarCorr) and the BLUPs (ranef).
  ##   mle_fit    a coxph with a `focus:group` interaction: supplies the raw per-group MLE
  ##              and SE (the unshrunk analog of the random slope).
  ##   data       modelling data, for the informative-subtype filter (events on both sides
  ##              of `focus`).
  ## MAP point = pooled + BLUP; MAP CI from the empirical-Bayes posterior variance
  ## tau^2 (1 - s_k), s_k = tau^2 / (tau^2 + se_k^2). Subtypes with no usable MLE (aliased,
  ## se 0, or non-finite) get v = Inf, so s_k = 0 and the posterior falls back to the prior.
  ## Returns one row per informative group, sorted by descending MAP, with columns subtype,
  ## mle_est/lo/hi, map_est/lo/hi (all natural-log HR) and attributes `pooled` and `tau`.
  mu   <- as.numeric(fixef(coxme_fit)[focus])
  tau2 <- as.numeric(unlist(VarCorr(coxme_fit)))[1]
  re   <- ranef(coxme_fit)[[group]]
  blup <- re[, 1]; names(blup) <- rownames(re)

  co  <- summary(mle_fit)$coefficients
  pat <- paste0("^", focus, ":", group)
  co  <- co[grepl(pat, rownames(co)), , drop = FALSE]
  rownames(co) <- sub(pat, "", rownames(co))
  beta <- co[, "coef"]; se <- co[, "se(coef)"]

  ev <- data[[event]] == 1; pos <- data[[focus]] == 1; g <- data[[group]]
  lvls <- levels(droplevels(as.factor(g)))
  info <- lvls[(tapply(ev & pos, g, any)[lvls] %in% TRUE) &
               (tapply(ev & !pos, g, any)[lvls] %in% TRUE)]

  tab <- data.frame(subtype = info, stringsAsFactors = FALSE)
  tab$beta <- beta[tab$subtype]; tab$se <- se[tab$subtype]; tab$blup <- blup[tab$subtype]
  tab$v <- ifelse(is.na(tab$beta) | !is.finite(tab$se) | tab$se == 0, Inf, tab$se^2)  # no MLE info -> prior
  tab$s <- tau2 / (tau2 + tab$v)
  map_sd <- sqrt(tau2 * (1 - tab$s))
  tab$mle_est <- tab$beta; tab$mle_lo <- tab$beta - 1.96 * tab$se; tab$mle_hi <- tab$beta + 1.96 * tab$se
  tab$map_est <- mu + tab$blup
  tab$map_lo  <- tab$map_est - 1.96 * map_sd; tab$map_hi <- tab$map_est + 1.96 * map_sd
  chk <- abs(tab$blup - tab$s * (tab$beta - mu))               # coxme BLUP vs EB-implied shrinkage
  tab <- tab[order(-tab$map_est), ]

  cat(sprintf("%s: pooled HR = %.2f; tau = %.3f; %d informative subtypes; max|BLUP - EB| = %.3f\n",
              focus, exp(mu), sqrt(tau2), nrow(tab), max(chk[is.finite(chk)])))
  out <- tab[, c("subtype", "mle_est", "mle_lo", "mle_hi", "map_est", "map_lo", "map_hi")]
  attr(out, "pooled") <- mu; attr(out, "tau") <- sqrt(tau2)
  out
}

caterpillar_forest <- function(df, pooled = attr(df, "pooled"),
                               xlim = NULL, ticks_at = NULL, span = 3.5,
                               base_size = 7, ci_width = 50,
                               col_mle = "#762a83", col_map = "#1b7837",
                               unstable_logwidth = 3) {
  ## Paired per-subtype caterpillar of an exposure log-HR: raw coxph MLE (interaction
  ## model) vs shrunk coxme MAP (BLUP), overlaid on one shared log-HR axis.
  ##
  ## Parameters
  ##   df      pre-sorted data frame, one row per subtype (rows are drawn top-to-bottom
  ##           in the order given), with columns:
  ##             subtype                 - subtype label (character)
  ##             mle_est, mle_lo, mle_hi - raw MLE point and 95% CI, natural-log-HR
  ##                                       scale; NA where not estimable
  ##             map_est, map_lo, map_hi - shrunk MAP point and 95% CI, log-HR scale
  ##   pooled  pooled mean exposure effect (log-HR; the coxme fixed effect), drawn as the
  ##           dashed vertical "shrinkage target" line (the reference line of the plot).
  ##           Defaults to attr(df, "pooled"), which caterpillar_table() sets.
  ##   xlim    plot window on the HR scale; estimates outside it are clamped to an
  ##           arrow at the edge so a single wild MLE cannot blow up the axis. Default
  ##           NULL builds a window geometrically centred on `pooled` (so the dashed
  ##           reference line sits at the horizontal centre rather than at HR = 1).
  ##   ticks_at  HR positions for the x-axis ticks (axis is log-spaced). Default NULL
  ##           uses the "nice" values in c(0.25,0.5,1,2,4,8) that fall inside xlim.
  ##   span    multiplicative half-width of the default window, xlim = exp(pooled) *
  ##           c(1/span, span). Ignored when xlim is supplied explicitly.
  ##   base_size base font size in points.
  ##   ci_width  width, in space characters, of the blank column that hosts the forest
  ##             CIs; larger = wider plotting panel (mirrors forest_with_strata's ci_width).
  ##   col_mle, col_map  colours for the MLE and MAP series.
  ##   unstable_logwidth  if a finite MLE's CI spans more than this many log-HR units
  ##             (an exp(unstable_logwidth)-fold HR range, ~20x at the default of 3),
  ##             the fit is treated as numerically unstable from (near-)perfect
  ##             separation: its text is shown as "unstable" while the point/CI still
  ##             draws as a clamped arrow.
  hr <- function(x) exp(x)
  if (is.null(pooled)) stop("caterpillar_forest(): supply `pooled`, or a `pooled` attribute on df (caterpillar_table sets it)")
  if (is.null(xlim)) xlim <- hr(pooled) * c(1 / span, span)   # geometrically centred on pooled
  if (is.null(ticks_at)) {
    cand <- c(0.25, 0.5, 1, 2, 4, 8)
    ticks_at <- cand[cand >= xlim[1] * 0.999 & cand <= xlim[2] * 1.001]
  }
  lo_cap <- xlim[1] / 3; hi_cap <- xlim[2] * 3
  clamp <- function(x) pmin(pmax(hr(x), lo_cap), hi_cap)   # keep >0 & near window so arrows draw

  fmt_ci <- function(e, l, h) {
    out <- sprintf("%.2f (%.2f–%.2f)", hr(e), hr(l), hr(h))
    unstable <- is.finite(e) & is.finite(l) & is.finite(h) & ((h - l) > unstable_logwidth)
    out[unstable] <- "unstable"
    out[is.na(e)] <- "not estimable"
    out
  }

  d <- df
  d$`Molecular subtype` <- d$subtype
  d$` ` <- strrep(" ", ci_width)
  d$`MLE HR (95% CI)` <- fmt_ci(d$mle_est, d$mle_lo, d$mle_hi)
  d$`MAP HR (95% CI)` <- fmt_ci(d$map_est, d$map_lo, d$map_hi)
  tab <- d[, c("Molecular subtype", " ", "MLE HR (95% CI)", "MAP HR (95% CI)")]

  tm <- forest_theme(base_size = base_size,
                     ci_pch = c(15, 16), ci_col = c(col_mle, col_map), ci_lwd = 1.4,
                     legend_name = "Estimate",
                     legend_value = c("MLE", "MAP"),
                     refline_gp = gpar(col = "grey55", lty = 2))

  forest(tab,
         est   = list(clamp(d$mle_est), clamp(d$map_est)),
         lower = list(clamp(d$mle_lo),  clamp(d$map_lo)),
         upper = list(clamp(d$mle_hi),  clamp(d$map_hi)),
         ci_column = 2, ref_line = hr(pooled),
         xlim = xlim, ticks_at = ticks_at, x_trans = "log", nudge_y = 0.25,
         theme = tm)
}

## named colour palette for KM groups, keyed on the group label (order-independent).
km_palette <- c(
  'ecDNA-'           = 'blue',
  'ecDNA+'           = 'red',
  'no amplification' = 'dodgerblue',
  'chromosomal'      = 'magenta',
  'ecDNA'            = 'red',
  'chr+ MYC-'        = 'orchid1',
  'chr+ MYC+'        = 'orchid4',
  'ecDNA+ MYC-'      = 'indianred1',
  'ecDNA+ MYC+'      = 'red4'
)

km_plot <- function(survObj, palette = km_palette){
  ## perform a KM analysis and generate the plot.
  ## Colours are looked up from `palette` by group label (the strata value with any "var="
  ## prefix stripped), in the actual strata order, so they track the group regardless of
  ## how the underlying factor is releveled.
  grp <- str_replace(names(survObj$strata), "^.*=", "")   # group labels, in strata order
  missing_grp <- setdiff(grp, names(palette))
  if (length(missing_grp))
    stop("km_plot: no colour defined for group(s): ", paste(missing_grp, collapse = ", "),
         ". Add them to `palette`.")
  colors <- unname(palette[grp])

  plt <- survObj %>%
   ggsurvfit(linewidth=0.5) +
   labs(x = 'Follow-up time (Months)',
        y = 'Overall Survival') +
   scale_color_manual(values = colors,
                      labels = grp) +
   scale_fill_manual(values = colors,
                     labels = grp) +
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
  return(plt)
}

save_ggplot <- function(outfile, width = 3, height = 3.5){
    ## Save the last displayed ggplot (ggsave's default, i.e. last_plot()) to out/ as png + svg.
    ## Defaults match the KM figure size; pass width/height for wider forests etc.
    pngName = paste(outfile, ".png", sep="")
    svgName = paste(outfile, ".svg", sep = "")
    ggsave(path="out", device="png", filename=pngName, width=width, height=height, units='in')
    ggsave(path="out", device="svg", filename=svgName, width=width, height=height, units='in')
}