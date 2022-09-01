###---------------------------------------------
# ME.440.825 Quantitative Neurogenomics
# Problem Set 1
# Dependencies and supplementary functions
# Author: Alina Spiegel
# Date: 8/30/2021
###---------------------------------------------

# load dependencies
library(ggplot2)

# function to calculate the statistical mode
# arguments: 
#   x         a vector or matrix 
getMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# function to create side-by-side box plots of gene expression for each
# cell class
# arguments:
#   df        a dataframe containing the column "cellClass" and a column 
#             corresponding to the argument geneName with values 
#             log2(expression + 1) for that gene
#   geneName  a string specifying the gene of interest. Must also be a column in
#             the dataframe df
makeBP <- function(df, geneName) {
  p<-ggplot(df, aes_string("cellClass", geneName)) + 
    geom_boxplot() +
    ylab(paste(geneName, 'log2(counts+1)'))
  
  print(p)
}
