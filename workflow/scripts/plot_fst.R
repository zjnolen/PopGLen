sink(file(snakemake@log[[1]], open="wt"), type = "message")

plot_fst <- function(fst_table, plotout) {
	require(ggplot2)
	require(ggtext)
	require(RColorBrewer)

	fst <- as.data.frame(read.table(fst_table, header = TRUE))
	fst2 <- fst
	fst2$pop1 <- fst$pop2
	fst2$pop2 <- fst$pop1
	fst <- rbind(fst,fst2)

	ggplot(fst, aes(pop1, pop2, fill=weight.fst, label=round(weight.fst,4))) +
		geom_tile() +
		geom_text() +
		labs(fill = "Weighted *F*~ST~") +
		scale_fill_distiller(palette = "YlOrRd", direction = 1) +
		theme_classic() +
		theme(
			aspect.ratio = 1,
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
			legend.title = element_markdown()
		)

	npop <- length(unique(fst$pop1))

	if (npop > 6) {
		wh <- 7+(0.77*npop)
	} else {
		wh <- 7
	}

	ggsave(plotout, width = wh, height = wh)

}

plot_fst(snakemake@input[[1]],snakemake@output[[1]])