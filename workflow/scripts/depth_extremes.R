# calculates upper and lower percentiles of lower depth, excluding sites 
# with depth = 0

sink(file(snakemake@log[[1]], open="wt"), type = "message")

chrom.dp.table <- read.table(snakemake@input[[1]], header = FALSE)
genome.dp <- colSums(chrom.dp.table)
genome.dp[1] <- 0
df <- data.frame(matrix(nrow = length(genome.dp), ncol = 0))
df$dp <- seq(from = 0, to = length(genome.dp)-1)
df$count <- genome.dp
quants <- with(df, quantile(rep(x = dp, times = count),
			probs = c(snakemake@params[["lower"]],
			snakemake@params[["upper"]])))
quants <- unname(quants)
mean <- with(df, mean(rep(x = dp, times = count)))
toprint <- c(mean,as.integer(quants))
write(toprint, snakemake@output[[1]])
