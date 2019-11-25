library(phyloseq)
library(dplyr)
load('Floresphy.RData')

phy_taxfilt <- filter(phy, 'body_site', "UBERON:feces")
save(phy_taxfilt, file='Floresphyfilt.RData')
