library(phyloseq)
library(dplyr)
load('Bolivarphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "feces")
save(phy_taxfilt, file='Bolivarphyfilt.RData')
