library(phyloseq)
library(dplyr)
load('GGMPphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "feces")
save(phy_taxfilt, file='GGMPphyfilt.RData')
