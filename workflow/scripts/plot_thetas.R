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


  # Some variables to automatically adjust the width of the plot based on the
  # number of categories on the x-axis. This is kind of based on the assumption
  # that 5ish fit comfortably at the default dimensions (7x7in), so we add 1/5
  # of the width for every sample above 5.
  npops <- length(unique(pestPGcomb$pop))
  if (npops <= 5) {
    wpop <- 7
  } else {
    wpop <- (7/5)*npops
  }

  ggplot(pestPGcomb, aes(x=pop, y=watterson)) +
    geom_violin(fill = "darkgoldenrod2") +
    stat_summary(fun.data = mean_cl_boot, geom = "pointrange", 
      color = "black") +
    labs(y = "Watterson's \316\270", x = "Population") +
    theme_classic()
  
  ggsave(paste0(plotpre,".watterson.pdf"), width = wpop, units = "in")

  ggplot(pestPGcomb, aes(x=pop, y=pi)) +
    geom_violin(fill = "darkgoldenrod2") +
    stat_summary(fun.data = mean_cl_boot, geom = "pointrange", 
      color = "black") +
    labs(y = "Nucleotide Diversity (\317\200)", x = "Population") +
    theme_classic()
  
  ggsave(paste0(plotpre,".pi.pdf"), width = wpop, units = "in")

  ggplot(pestPGcomb, aes(x=pop, y=Tajima)) +
    geom_violin(fill = "darkgoldenrod2") +
    geom_hline(yintercept=0, linetype = "dashed") +
    stat_summary(fun.data = mean_cl_boot, geom = "pointrange", 
      color = "black") +
    labs(y = "Tajima's D", x = "Population") +
    ylim(c(-2,2)) +
    theme_classic()
  
  ggsave(paste0(plotpre,".tajima.pdf"), width = wpop, units = "in")

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