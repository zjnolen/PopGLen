sink(file(snakemake@log[[1]], open="wt"), type = "message")

source("/usr/local/bin/visFuns.R")

plot_evaladmix <- function(pop,qopt,corres,k) {

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
		title=paste0("Evaluation of admixture proportions with K=",k),
		max_z=0.1, min_z=-0.1)
	
	return(plot)

}

pdf(snakemake@output[[1]], width = 11.8, height = 9.8)
plot_evaladmix(snakemake@input[["pops"]],snakemake@input[["qopt"]],
	snakemake@input[["corres"]],snakemake@wildcards[["kvalue"]])
