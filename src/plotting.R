Sys.setenv(LANGUAGE = "en")

library(ggplot2)
library(svglite)
library(ragg)

save_ggplot <- function(outfile, plot = last_plot(), path="out", width = 3, height = 3.5, units="in"){
    ggsave(path=path, filename=paste0(outfile, ".png"), plot=plot,
           device=ragg::agg_png,   width=width, height=height, units=units)
    ggsave(path=path, filename=paste0(outfile, ".svg"), plot=plot,
           device=svglite::svglite, width=width, height=height, units=units)
}

save_base_plot <- function(outfile, plot_fun, path="out", width = 3, height = 3.5, pointsize = 7, family = "Arial"){
    ## Save a base R plot as png and svg, matching save_ggplot and the Arial 7 pt figure style.
    ## plot_fun: a function with no arguments that draws the plot, eg. function() plotCorrMatrix(C). It is called once per device.
    ## Text is drawn at pointsize; functions that scale their own text (eg. heatmap.2 cexRow/cexCol) need those set accordingly.
    ## width, height: inches.
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    devices <- list(
        png = function(file) ragg::agg_png(file, width = width, height = height, units = "in", res = 300, pointsize = pointsize),
        svg = function(file) svglite::svglite(file, width = width, height = height, pointsize = pointsize))
    for (ext in names(devices)) {
        devices[[ext]](file.path(path, paste0(outfile, ".", ext)))
        par(family = family)
        plot_fun()
        dev.off()
    }
}