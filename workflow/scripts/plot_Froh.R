sink(file(snakemake@log[[1]], open="wt"), type = "message")

aggregate_roh <- function(rohlist,poplist,samplelist,lenautos) {
  samples <- as.data.frame(read.table(samplelist, header=TRUE))
  Froh <- data.frame()
  
  for (i in 1:length(poplist)) {

    for (l in c(100000,1000000,2500000)) {
      if (file.size(rohlist[i]) == 0) {
        roh <- data.frame(matrix(ncol = 5, nrow = 0))
      } else {
        roh <- as.data.frame(read.table(rohlist[i], header = FALSE))
      }
      for (s in samples$sample[samples$population==poplist[i]]) {
        sampleroh <- (subset(roh, roh[,4] == s))
        nroh <- nrow(sampleroh[sampleroh[,5]>l,])
        sumroh <- sum(sampleroh[,5][sampleroh[,5] > l])
        Fsample <- sumroh/lenautos
        row <- c(s,sumroh,Fsample,nroh,poplist[i],l)
        Froh <- rbind(Froh, row)
      }

    }

  }

  colnames(Froh) <- c("sample","lenroh","Froh","Nroh","population","minroh")

  Froh$sample <- as.factor(Froh$sample)
  Froh$population <- as.factor(Froh$population)
  Froh$lenroh <- as.integer(Froh$lenroh)
  Froh$Froh <- as.numeric(Froh$Froh)
  Froh$Nroh <- as.integer(Froh$Nroh)
  Froh$minroh <- as.factor(Froh$minroh)

  return(Froh)
}

plot_roh <- function(aggroh,plotpre) {
  require(ggplot2)

  ggplot(aggroh, aes(x=population,y=Froh)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    scale_y_continuous(limits = c(0, NA)) +
    facet_wrap(~ minroh) +
    theme_classic()
  
  ggsave(paste0(plotpre,".froh.pdf"))
  
  # ggplot(aggroh, aes(x=Nroh,y=lenroh)) +
  #   geom_point(aes(col=population)) +
  #   geom_smooth(method='lm') +
  #   theme_classic() +
  #   facet_grid(cols = vars(minroh))

  # ggsave(paste0(plotpre,".rohreg.pdf"))

  write.table(aggroh, file = paste0(plotpre,".froh.tsv"), quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = FALSE)

}

aggroh <- aggregate_roh(snakemake@input[["roh"]],
  snakemake@params[["popnames"]],
  snakemake@input[["inds"]],
  read.table(snakemake@input[["autos"]], sep = "\t")[,2])

plot_roh(aggroh,snakemake@params[["outpre"]])