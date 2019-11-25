library(phyloseq)
library(dplyr)
load('TZphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='TZphyfilt.RData')
