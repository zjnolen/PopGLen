# Reads a 2D SFS of two individuals and estimates R0, R1, and KING-robust
# kinship coefficients as in Waples et al. 2019, MolEcol. R code adapted from 
# supplementary materials available at:
# https://github.com/rwaples/freqfree_suppl/blob/master/read_realSFS.R

calc_coeffs <- function(filepath){
  df = read.csv(filepath, sep =' ', header = FALSE)
  colnames(df) = c('A', 'D', 'G', 'B', 'E', 'H', 'C', 'F', 'I', 'drop')
  df['drop'] = NULL
  
  df['HETHET'] = df['E']
  df['IBS0'] = df['C'] + df['G']
  df['IBS1'] = df['B'] + df['D'] + df['F'] + df['H']
  df['IBS2'] = df['A'] + df['E'] + df['I']
  
  coeffs = data.frame(matrix(nrow=1,ncol=0))
  # the derived stats 
  coeffs['R0'] = df['IBS0'] / df['HETHET']
  coeffs['R1'] = df['HETHET'] / (df['IBS0'] +  df['IBS1'])
  # KING-robust kinship
  coeffs['Kin'] = 
    (df['HETHET'] - 2*(df['IBS0'])) / (df['IBS1'] + 2*df['HETHET'])
  return(round(coeffs, digits = 8))
}

toprint = data.frame(matrix(nrow=1,ncol=0))
toprint['ind1'] = snakemake@wildcards[["ind1"]]
toprint['ind2'] = snakemake@wildcards[["ind2"]]
toprint = cbind(toprint,calc_coeffs(snakemake@input[[1]]))
write.table(toprint, file = snakemake@output[[1]], sep = "\t", 
  row.names = FALSE, col.names = FALSE, quote = FALSE)