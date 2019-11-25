library(phyloseq)
library(dplyr)
load('Infantphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', 'Stool')
save(phy_taxfilt, file='Infantphyfilt.RData')
