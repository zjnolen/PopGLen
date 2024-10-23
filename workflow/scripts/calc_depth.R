sink(file(snakemake@log[[1]], open="wt"), type = "message")

total_len <- function(filtersum) {
	df <- read.table(filtersum, header = FALSE, sep = "\t")
	total <- df[nrow(df), 2]
	return(total)
}

summarize_depth <- function(depth, sample, length) {
	depth <- scan(depth)
	depth[1] <- 0
	depth[1] <- length-sum(depth)
	df <- data.frame(matrix(nrow=length(depth), ncol = 0))
	df$dp <- seq(from = 0, to = length(depth)-1)
	df$count <- depth

	mean <- round(sum(df$dp * df$count)/sum(df$count), digits = 4)
	stdev <- round(
		sqrt(sum(df$count*(df$dp-mean)^2)/(sum(df$count)-1)),
		digits = 4
	)

	return(c(sample,mean,stdev))
}



write(
	summarize_depth(
		snakemake@input[["sample_hist"]], 
		snakemake@wildcards[["sample"]],
		total_len(snakemake@input[["filtersum"]])
	),
	file = snakemake@output[["sample_summ"]],
	sep = "\t",
	ncolumns=3
)