sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(dplyr)
library(ggplot2)


ibstable_raw <- read.table(snakemake@input[["ibs"]], sep = "\t", header = TRUE)
poplist <- read.table(snakemake@input[["pops"]], sep = "\t", header = TRUE)
plotpre <- snakemake@params[["plotpre"]]

ibstable <- merge(ibstable_raw, poplist, by = "sample")

ggplot(ibstable, aes(x = population, y = ibs.to.ref)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  ylab("Proportion sites identical to reference") +
  xlab("Population") +
  theme_classic()

ggsave(paste0(plotpre, ".population.svg"))

ggplot(ibstable, aes(x = time, y = ibs.to.ref)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  ylab("Proportion sites identical to reference") +
  xlab("Time Period") +
  theme_classic()

ggsave(paste0(plotpre, ".time.svg"))

ggplot(ibstable, aes(x = depth, y = ibs.to.ref)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  ylab("Proportion sites identical to reference") +
  xlab("Depth Class") +
  theme_classic()

ggsave(paste0(plotpre, ".depth.svg"))

write.table(ibstable, file = snakemake@output[["table"]], quote=FALSE,
  sep = '\t', row.names = FALSE)