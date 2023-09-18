# calculates upper and lower percentiles of lower depth
# adapted from https://github.com/KHanghoj/leopardpaper/blob/master/mapping_and_qc/reference_filter/depth/04_quantile_thresholds.R

sink(file(snakemake@log[[1]], open="wt"), type = "message")

chrom.dp.table <- read.table(snakemake@input[[1]], header = FALSE)
genome.dp <- colSums(chrom.dp.table)
genome.dp.cumsum <- cumsum(as.numeric(genome.dp))

#qup <- genome.dp.cumsum[length(genome.dp.cumsum)]*snakemake@params[["upper"]]
#qlow <- genome.dp.cumsum[length(genome.dp.cumsum)]*snakemake@params[["lower"]]
#upper <- min(which(genome.dp.cumsum > qup))
#lower <- min(which(genome.dp.cumsum > qlow))

df <- data.frame(matrix(nrow = length(genome.dp), ncol = 0))
df$dp <- seq(from = 0, to = length(genome.dp)-1)
df$count <- genome.dp
mean <- with(df, mean(rep(x = dp, times = count)))
median <- with(df, median(rep(x = dp, times = count)))
quants <- c(
  median * snakemake@params[["lower"]],
  median * snakemake@params[["upper"]]
)
toprint <- c(mean,as.integer(quants),median)
write(toprint, snakemake@output[[1]])