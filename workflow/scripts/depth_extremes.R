# calculates upper and lower percentiles of lower depth
# adapted from https://github.com/KHanghoj/leopardpaper/blob/master/mapping_and_qc/reference_filter/depth/04_quantile_thresholds.R

sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(ggplot2)

chrom.dp.table <- read.table(snakemake@input[[1]], header = FALSE)
genome.dp <- colSums(chrom.dp.table)
df <- data.frame(matrix(nrow = length(genome.dp), ncol = 0))
df$dp <- seq(from = 0, to = length(genome.dp)-1)
df$count <- genome.dp
mean <- sum(df$dp * df$count)/sum(df$count)
genome.dp.cumsum <- cumsum(as.numeric(genome.dp))
qmed <- genome.dp.cumsum[length(genome.dp.cumsum)]*0.5
median <- min(which(genome.dp.cumsum >= qmed))-1

if (snakemake@params[["method"]] == "percentile") {
  qup <- genome.dp.cumsum[length(genome.dp.cumsum)]*snakemake@params[["upper"]]
  qlow <- genome.dp.cumsum[length(genome.dp.cumsum)]*snakemake@params[["lower"]]
  upper <- min(which(genome.dp.cumsum > qup))-1
  if (snakemake@params[["lower"]] == 0) {
    lower <- 0
  } else {
    lower <- max(min(which(genome.dp.cumsum > qlow))-1, 0)
  }
  quants <- c(lower, upper)
} else if (snakemake@params[["method"]] == "median") {
  quants <- c(
    median * snakemake@params[["lower"]],
    median * snakemake@params[["upper"]]
  )
}

## Make plot of depth distribution and cutoffs
xlim <- min(
  which(genome.dp.cumsum > genome.dp.cumsum[length(genome.dp.cumsum)]*0.995)
)
ggplot(df[c(1:xlim), ], aes(y = count, x = dp)) +
  geom_bar(stat='identity', fill = "gray40") +
  ggtitle(paste0(
    "Global depth distribution ",
    snakemake@wildcards[["population"]],
    " depth samples"
    )
  ) +
  xlab("Global Depth") +
  ylab("Count") +
  xlim(c(0,xlim)) +
  geom_vline(xintercept=quants[1], color = "red") +
  geom_vline(xintercept=quants[2], color = "red") +
  theme_classic()

ggsave(snakemake@output[["plot"]])

toprint <- c(mean,as.integer(quants),median)
write(toprint, snakemake@output[["summ"]])