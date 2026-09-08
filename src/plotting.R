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
