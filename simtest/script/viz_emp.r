#!/usr/bin/env Rscript
library(ggplot2)
library(readr)
library(dplyr)

# theme_set(theme_bw())
# theme_set(theme_classic())
# theme_set(theme_light()())
options("width"=200)

dir.create("../viz")

results <- read_csv("results_empirical.csv")

# get overall statistics

round(mean(results$abs_err_mean), digits=2)
formatC(mean(results$norm_err_mean), format = "e", digits = 1)

diffs <- aggregate(results[, c("abs_err_mean", "norm_err_mean")], list(rarify_fract=results$rarify_fract, query_fract=results$query_fract), mean)

diffs

# data.frame(diff(as.matrix(diffs)))

diffs <- mutate(diffs,
  abs_err_diff=((abs_err_mean / min(abs_err_mean))-1)*100,
  norm_err_diff=((norm_err_mean / min(norm_err_mean))-1)*100
)

diffs

diffs <- aggregate(results[, c("abs_err_mean", "norm_err_mean")], list(rarify_fract=results$rarify_fract), mean)

# diffs

diffs <- mutate(diffs,
  abs_err_diff=((abs_err_mean / min(abs_err_mean))-1)*100,
  norm_err_diff=((norm_err_mean / min(norm_err_mean))-1)*100
)

diffs

diffs <- aggregate(results[, c("abs_err_mean", "norm_err_mean")], list(query_fract=results$query_fract), mean)

# diffs

diffs <- mutate(diffs,
  abs_err_diff=((abs_err_mean / min(abs_err_mean))-1)*100,
  norm_err_diff=((norm_err_mean / min(norm_err_mean))-1)*100
)

diffs

params <- list("rarify_fract")
# metrics <- list("rel_err_mean", "rel_norm_err_mean", "norm_norm_krd")
metrics <- list("count_diff")
modes <- list("rootings","bootstrap","outgroup")

for (param in params) {
  for (metric in metrics) {
    p <- ggplot( results, aes_string(x=param) )
    for (mode in modes) {
      slice=results[results$scrapp_mode == mode,]
      p <- p + stat_summary( data=slice,
          aes_string(y=metric, group="1", linetype="scrapp_mode", colour="scrapp_mode"),
          fun.y=mean, geom="line",group=1 ) +
      stat_summary( data=slice,
          aes_string(y=metric, group="1", colour="scrapp_mode"),
          fun.y = mean, geom = "point", group=1, show.legend = FALSE) +
      stat_summary( data=slice,
          aes_string(y=metric, group="1", linetype="scrapp_mode", colour="scrapp_mode"),
          fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", group=1, show.legend = FALSE)
      # geom_errorbar(yintercept = "mean",width=0.8,aes(ymax=..metric..,ymin=..metric..))
    }
    ggsave(paste("../viz/",param,"|",metric,".pdf", sep=""), plot=p, device="pdf")
  }
}
