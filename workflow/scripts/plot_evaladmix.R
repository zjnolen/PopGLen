sink(file(snakemake@log[[1]], open="wt"), type = "message")

plot_evaladmix <- function(pop,qopt,corres,k) {

	source("https://raw.githubusercontent.com/GenisGE/evalAdmix/2a51aebaca70c9d3b5fb359d6c5c40145c58fce5/visFuns.R")

	# read population labels and estimated admixture proportions
	pop<-read.table(pop, header = TRUE, sep = "\t")
	pop
	q<-read.table(qopt)
	q

	# order according to population
	ord<-orderInds(pop = as.vector(pop[,2]), q = q)

	# read in correlation matrix
	r<-as.matrix(read.table(corres))

	# Plot correlation of residuals
	plot <- plotCorRes(cor_mat = r, pop = as.vector(pop[,2]), ord=ord, 
		title=paste0("Evaluation of admixture proportions with K=",k))
	
	return(plot)

}

png(snakemake@output[[1]], res = 300, width = 30, height = 25, units = "cm")
plot_evaladmix(snakemake@input[["pops"]],snakemake@input[["qopt"]],
	snakemake@input[["corres"]],snakemake@wildcards[["kvalue"]])
dev.off()