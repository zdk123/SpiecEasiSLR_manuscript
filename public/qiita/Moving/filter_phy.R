library(phyloseq)
library(dplyr)
load('Movingphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='Movingphyfilt.RData')
