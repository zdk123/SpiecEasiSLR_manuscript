library(phyloseq)
library(dplyr)
load('DHMphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "Stool")
save(phy_taxfilt, file='DHMphyfilt.RData')
