library(phyloseq)
library(dplyr)
load('AGPphy.RData')

phy_taxfilt <- filter(phy)
save(phy_taxfilt, file='AGPphyfilt.RData')
