# calculates upper and lower percentiles of lower depth, excluding sites 
# with depth = 0

chrom.dp.table <- read.table(snakemake@input[[1]], header = FALSE)
genome.dp <- colSums(chrom.dp.table)
genome.dp[1] <- 0
df <- data.frame(matrix(nrow = length(genome.dp), ncol = 0))
df$dp <- seq(from = 0, to = length(genome.dp)-1)
df$count <- genome.dp
quants <- with(df, quantile(rep(x = dp, times = count),
			probs = c(snakemake@params[["lower"]],
			snakemake@params[["upper"]])))
quants
quants <- unname(quants)
mean <- with(df, mean(rep(x = dp, times = count)))
mean
toprint <- c(mean,as.integer(quants))
toprint
write(toprint, snakemake@output[[1]])
