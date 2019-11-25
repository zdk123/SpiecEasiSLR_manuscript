library(phyloseq)
library(dplyr)
load('Hadzaphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "feces")
save(phy_taxfilt, file='Hadzaphyfilt.RData')
