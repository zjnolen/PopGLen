sink(file(snakemake@log[[1]], open="wt"), type = "message")

plot_fst <- function(fst_table, plotout) {
	require(ggplot2)

	fst <- as.data.frame(read.table(fst_table, header = TRUE))
	fst2 <- fst
	fst2$pop1 <- fst$pop2
	fst2$pop2 <- fst$pop1
	fst <- rbind(fst,fst2)

	ggplot(fst, aes(pop1,pop2,fill=weight.fst,label=round(weight.fst,4))) +
		geom_tile() +
		geom_text() +
		theme_classic() +
		theme(aspect.ratio = 1)

	ggsave(plotout)

}

plot_fst(snakemake@input[[1]],snakemake@output[[1]])