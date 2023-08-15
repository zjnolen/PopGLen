sink(file(snakemake@log[[1]], open="wt"), type = "message")

aggregate_heterozygosity <- function(sfslist,pop) {

  s <- as.data.frame(read.table(pop, header = TRUE))

  hz <- data.frame(matrix(NA, ncol = 3, nrow = 0))
  colnames(hz) <- c("sample","pop","heterozygosity")

  for (i in 1:length(sfslist)) {
    
    sfs <- scan(sfslist[i])
    samplehz <- sfs[2]/sum(sfs)
    hz[nrow(hz) + 1,] <- c(s$sample[i],s$population[i],samplehz*1000)
  }

  hz$pop <- as.factor(hz$pop)
  hz$pop <- factor(hz$pop, levels = sort(unique(hz$pop)))
  hz$heterozygosity <- as.numeric(hz$heterozygosity)
  return(hz)

}

plot_heterozygosity <- function(agghz, plotout) {
  
  require(ggplot2)

  ggplot(agghz, aes(x=pop, y=heterozygosity, fill=pop)) +
    geom_boxplot() +
    xlab("Population") +
    ylab("Heterozygous sites per 1000bp") +
    theme_classic()

  ggsave(plotout)

}

agghz <- aggregate_heterozygosity(
  snakemake@input[["sfs"]],snakemake@input[["popfile"]])

write.table(agghz, file = snakemake@output[[1]], row.names=FALSE, sep = "\t",
              quote=FALSE)

plot_heterozygosity(agghz,snakemake@output[[2]])