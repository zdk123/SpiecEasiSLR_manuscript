library(phyloseq)
library(dplyr)
load('ECAMphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "Stool|Rectal")
save(phy_taxfilt, file='ECAMphyfilt.RData')
