library(phyloseq)
library(dplyr)
load('HMPphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='HMPphyfilt.RData')
