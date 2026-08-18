sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(data.table)
library(dplyr)

# Reads in an LD table from ngsLD and prunes variants for which the maximum
# r2 between upstream SNPs within a give distance exceeds a threshold.
# Mimics behavior of bcftools +prune -w <window_size> -m <r2_thresh>, but for
# ngsLD outputs. ngsLD must have been run with the max_kb_dist set to the size
# of the window you want to use. Output is a list of sites in approximate
# linkage equilibrium, but they won't be in order.

ld_table <- snakemake@input[["ld"]]
pos_file <- snakemake@input[["pos"]]
out_list <- snakemake@output[["sites"]]
r2_thresh <- as.numeric(snakemake@wildcards[["r2"]])

if (file.size(ld_table) == 0) {

  file.create(out_list)

} else {

  ld <- fread(ld_table, sep = "\t", header = FALSE)

  maxr2s <- ld %>%
    group_by(V2) %>%
    summarize(maxr2 = max(V7))

  linked_snps <- maxr2s[maxr2s$maxr2 > r2_thresh, ]$V2

  pos <- fread(pos_file, sep = "\t", header = FALSE)

  pos <- paste(pos$V1, pos$V2, sep = ":")

  unlinked_snps <- pos[!pos %in% linked_snps]

  write.table(
    unlinked_snps,
    file = out_list,
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )

}
