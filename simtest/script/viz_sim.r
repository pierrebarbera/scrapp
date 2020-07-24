#!/usr/bin/env Rscript
library(ggplot2)
library(readr)
library(egg)

out.path <- "../viz"
if ( !dir.exists(out.path) ) {dir.create(out.path, recursive=TRUE)}

params <- list("species", "mut_rate", "pop_size", "prune_fract", "sample_size", "seq_length")
params_in_log <- list("mut_rate", "pop_size")
# params <- list("pop_size", "species", "prune_fract")
# params <- list("species")
# metrics <- list("norm_err_mean", "abs_err_mean", "rel_err_mean", "rel_norm_err_mean", "norm_norm_krd", "norm_krd", "krd")

metrics <- list("abs_err_mean", "rel_norm_err_mean", "null_fract")
mtolabel <- c("MAE", "NMRE", "% of results are null model")

# metrics <- list("null_fract")
# mtolabel <- c("null_fract")

names(mtolabel) <- metrics

modes <- list("rootings","bootstrap","outgroup")

prefix <- "new_test"

plot_list = list()

for (param in params) {
  results <- read_csv(file.path( "results", prefix, paste(param,".csv", sep="") ) )
  for (metric in metrics) {
    p <- ggplot( results, aes_string(x=param) ) +
      ylab(mtolabel[[ metric ]])
    for (mode in modes) {
      slice=results[ results$scrapp_mode == mode, ]
      p <- p + stat_summary( data=slice,
          aes_string(y=metric, group="1", linetype="scrapp_mode", colour="scrapp_mode"),
          fun.y=mean, geom="line",group=1) +
      stat_summary( data=slice,
          aes_string(y=metric, group="1", colour="scrapp_mode"),
          fun.y = mean, geom = "point", group=1, show.legend = FALSE) +
      stat_summary( data=slice,
          aes_string(y=metric, group="1", linetype="scrapp_mode", colour="scrapp_mode"),
          fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", group=1, show.legend = FALSE)
      p <- p + theme(legend.title=element_blank())
      if (param %in% params_in_log) {
        p <- p + scale_x_continuous(trans='log10')
      }
      # p <- p + scale_color_discrete(name="SCRAPP mode")

    plot_list.append( p )
    ggsave(paste("../viz/",param,"|",metric,".pdf", sep=""), plot=p, device="pdf", height=2.5, width=4)
  }
}


ggarrange(plot_list)

ggsave("../viz/complete.pdf", plot=plot_mat, device="pdf", height=4, width=4)


# width=5, height=4
