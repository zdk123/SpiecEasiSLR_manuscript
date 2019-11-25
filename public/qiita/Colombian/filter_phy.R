library(phyloseq)
library(dplyr)
load('Colombianphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "feces")
save(phy_taxfilt, file='Colombianphyfilt.RData')
