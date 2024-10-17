sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(dplyr)
library(ggplot2)


ibstable_raw <- read.table(snakemake@input[["ibs"]], sep = "\t", header = TRUE)
poplist <- read.table(snakemake@input[["pops"]], sep = "\t", header = TRUE)
plotpre <- snakemake@params[["plotpre"]]

# Some variables to automatically adjust the width of the plot based on the
# number of categories on the x-axis. This is kind of based on the assumption
# that 5ish fit comfortably at the default dimensions (7x7in), so we add 1/5 of
# the width for every sample above 5.
npops <- length(unique(poplist$population))
if (npops <= 5) {
  wpop <- 7
} else {
  wpop <- (7/5)*npops
}
ndp <- length(unique(poplist$depth))
if (ndp <= 5) {
  wdep <- 7
} else {
  wdep <- (7/5)*ndp
}

ibstable <- merge(ibstable_raw, poplist, by = "sample")

ggplot(ibstable, aes(x = population, y = ibs.to.ref)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  ylab("Proportion sites identical to reference") +
  xlab("Population") +
  theme_classic()

ggsave(paste0(plotpre, ".population.pdf"), width = wpop, units = "in")

ggplot(ibstable, aes(x = time, y = ibs.to.ref)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  ylab("Proportion sites identical to reference") +
  xlab("Time Period") +
  theme_classic()

ggsave(paste0(plotpre, ".time.pdf"))

ggplot(ibstable, aes(x = depth, y = ibs.to.ref)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  ylab("Proportion sites identical to reference") +
  xlab("Depth Class") +
  theme_classic()

ggsave(paste0(plotpre, ".depth.pdf"), width = wdep, units = "in")

write.table(ibstable, file = snakemake@output[["table"]], quote=FALSE,
  sep = '\t', row.names = FALSE)