# Reads IBS output from ANGSD and estimates R0, R1, and KING-robust
# kinship coefficients as in Waples et al. 2019, MolEcol. R code adapted from 
# supplementary materials available at:
# https://github.com/rwaples/freqfree_suppl/blob/master/read_IBS.R

sink(file(snakemake@log[[1]], open="wt"), type = "message")

library("stringr")

source(
  "https://raw.githubusercontent.com/rwaples/freqfree_suppl/a71b7a17fffa9ff50311540f0ec24db3b94e17a7/read_IBS.R"
)

inds <- read.table(snakemake@input[["inds"]], header = TRUE)

df <- c()

for (i in snakemake@input[["ibs"]]) {
  kins <- do_derived_stats(
    read_ibspair_model0(
      i
    )
  )
  df <- rbind(df, kins)
}

df <- aggregate(cbind(A,B,C,D,E,F,G,H,I) ~ pair, data = df, FUN = sum)
df['HETHET'] = df['E']
df['IBS0'] = df['C'] + df['G']
df['IBS1'] = df['B'] + df['D'] + df['F'] + df['H']
df['IBS2'] = df['A'] + df['E'] + df['I']
df['R0'] = df['IBS0'] / df['HETHET']
df['R1'] = df['HETHET'] / (df['IBS0'] +  df['IBS1'])
df['Kin'] = (df['HETHET'] - 2*(df['IBS0'])) / (df['IBS1'] + 2*df['HETHET'])

df[c('ind1','ind2')] <- str_split_fixed(df$pair, "_", 2)

df$ind1 <- as.numeric(df$ind1) + 1
df$ind2 <- as.numeric(df$ind2) + 1

df <- merge(df, inds, by.x = "ind1", by.y = 0)
df$ind1 <- df$sample
df$pop1 <- df$population
df <- df[,c("ind1","pop1","ind2","R0","R1","Kin")]
df <- merge(df, inds, by.x = "ind2", by.y = 0)
df$ind2 <- df$sample
df$pop2 <- df$population
df <- df[,c("ind1","ind2","pop1","pop2","R0","R1","Kin")]
df <- df[rev(order(df$Kin)),]


write.table(df, file = snakemake@output[[1]], sep = "\t", 
  row.names = FALSE, col.names = TRUE, quote = FALSE)