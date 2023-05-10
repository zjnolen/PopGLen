sink(file(snakemake@log[[1]], open="wt"), type = "message")

bed_total_len <- function(bedfile) {
	df <- data.frame(read.table(bedfile))
	total <- sum(df[,3]-df[,2])
	return(total)
}

summarize_depth <- function(depth, sample, length) {
	depth <- scan(depth)
	depth[1] <- 0
	depth[1] <- length-sum(depth)
	df <- data.frame(matrix(nrow=length(depth), ncol = 0))
	df$dp <- seq(from = 0, to = length(depth)-1)
	df$count <- depth

	mean <- round(with(df, mean(rep(x = dp, times = count))), digits = 8)
	stdev <- round(with(df, sd(rep(x=dp,times=count))), digits = 8)

	return(c(sample,mean,stdev))
}



write(
	summarize_depth(snakemake@input[["sample_hist"]], 
					snakemake@wildcards[["sample"]],
					bed_total_len(snakemake@input[["bed"]])
					),
	file = snakemake@output[["sample_summ"]],
	sep = "\t",
	ncolumns=3
)