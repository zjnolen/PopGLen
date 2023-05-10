sink(file(snakemake@log[[1]], open="wt"), type = "message")

combine_pestPG <- function(pestPGlist, popnames) {

  combined <- c()
  
  for (i in 1:length(pestPGlist)) {

    theta <- as.data.frame(read.table(pestPGlist[[i]], header = TRUE, 
      comment.char = ""))
    theta <- theta[!theta$nSites == 0,]
    theta$watterson <- theta$tW / theta$nSites
    theta$pi <- theta$tP / theta$nSites
    theta$pop <- popnames[i]
    combined <- rbind(combined,theta)

  }

  combined$pop <- as.factor(combined$pop)
  combined$pop <- factor(combined$pop, levels = sort(unique(combined$pop)))
  return(combined)

}

plot_thetas <- function(pestPGcomb, plotpre) {
  
  require(ggplot2)
  require(Hmisc)

  ggplot(pestPGcomb, aes(x=pop, y=watterson, fill=pop)) +
    geom_violin() +
    stat_summary(fun.data = mean_cl_boot, geom = "pointrange", 
      color = "black") +
    theme_classic()
  
  ggsave(paste0(plotpre,".watterson.pdf"))

  ggplot(pestPGcomb, aes(x=pop, y=pi, fill=pop)) +
    geom_violin() +
    stat_summary(fun.data = mean_cl_boot, geom = "pointrange", 
      color = "black") +
    theme_classic()
  
  ggsave(paste0(plotpre,".pi.pdf"))

  ggplot(pestPGcomb, aes(x=pop, y=Tajima, fill=pop)) +
    geom_violin() +
    stat_summary(fun.data = mean_cl_boot, geom = "pointrange", 
      color = "black") +
    ylim(c(-2,2)) +
    theme_classic()
  
  ggsave(paste0(plotpre,".tajima.pdf"))

}

combined <- combine_pestPG(snakemake@input,snakemake@params[["popnames"]])

plot_thetas(combined,snakemake@params[["outpre"]])