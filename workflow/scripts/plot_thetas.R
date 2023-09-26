sink(file(snakemake@log[[1]], open="wt"), type = "message")

combine_pestPG <- function(pestPGlist, popnames, minsites) {

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
  combined <- combined[combined$nSites >= minsites,]
  combined$Tajima <- sapply(
    combined$Tajima,
    function(x) replace(x, is.infinite(x), NA)
  )
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

theta_table <- function(pestPGcomb, tabpre) {

  require(Hmisc)
  
  pestPGcomb$watterson <- pestPGcomb$tW/pestPGcomb$nSites
  watterson <- aggregate(pestPGcomb$watterson, list(sites = pestPGcomb$pop), FUN=function(x) {smean.cl.boot(x, B = 2000, na.rm = TRUE)})
  watterson <- cbind(Population = watterson$sites, as.data.frame.matrix(watterson$x, row.names = NULL))
  write.table(watterson, file = paste0(tabpre,".watterson.mean.tsv"), quote=FALSE, sep = '\t', row.names = FALSE)

  pestPGcomb$pi <- pestPGcomb$tP/pestPGcomb$nSites
  pi <- aggregate(pestPGcomb$pi, list(sites = pestPGcomb$pop), FUN=function(x) {smean.cl.boot(x, B = 2000, na.rm = TRUE)})
  pi <- cbind(Population = pi$sites, as.data.frame.matrix(pi$x, row.names = NULL))
  write.table(pi, file = paste0(tabpre,".pi.mean.tsv"), quote=FALSE, sep = '\t', row.names = FALSE)


  tajima <- aggregate(pestPGcomb$Tajima, list(sites = pestPGcomb$pop), FUN=function(x) {smean.cl.boot(x, B = 2000, na.rm = TRUE)})
  tajima <- cbind(Population = tajima$sites, as.data.frame.matrix(tajima$x, row.names = NULL))
  write.table(tajima, file = paste0(tabpre,".tajima.mean.tsv"), quote=FALSE, sep = '\t', row.names = FALSE)

}

combined <- combine_pestPG(snakemake@input,snakemake@params[["popnames"]],snakemake@params[["minsites"]])

plot_thetas(combined,snakemake@params[["plotpre"]])
theta_table(combined,snakemake@params[["tabpre"]])