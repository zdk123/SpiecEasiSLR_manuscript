library(phyloseq)
library(dplyr)
load('Sandbornphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', 'stool')
save(phy_taxfilt, file='Sandbornphyfilt.RData')
