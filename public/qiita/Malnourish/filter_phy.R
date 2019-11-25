library(phyloseq)
library(dplyr)
load('Malnourishphy.RData')

phy_taxfilt <- filter(phy, 'body_site', "UBERON:feces")
save(phy_taxfilt, file='Malnourishphyfilt.RData')
