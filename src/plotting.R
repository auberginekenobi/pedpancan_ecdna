Sys.setenv(LANGUAGE = "en")

library(ggplot2)
library(svglite)
library(ragg)

COLORS_YML <- "../../config/colors.yml"  # relative to notebooks/<analysis>/

get_palette <- function(name, shade = NULL, path = COLORS_YML){
    ## Named character vector (label -> hex) for palette `name` in config/colors.yml,
    ## ready for scale_colour_manual(values = ...) / scale_fill_manual(values = ...).
    ## shade: for palettes with dark/light entries (tumor_type, tumor_subtype), which one to return.
    ##        Missing shades (YAML null) come back as NA.
    pal <- yaml::read_yaml(path)[[name]]
    if (is.null(pal)) stop("get_palette: no palette '", name, "' in ", path)
    if (!is.null(shade)) pal <- lapply(pal, function(x) x[[shade]])
    else if (any(vapply(pal, is.list, logical(1))))
        stop("get_palette: palette '", name, "' has shades; pass shade = 'dark' or 'light'")
    vapply(pal, function(x) if (is.null(x)) NA_character_ else x, character(1))
}

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