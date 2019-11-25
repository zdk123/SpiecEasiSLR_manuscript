library(phyloseq)
library(dplyr)
load('Babiesphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "Anal|Feces")
save(phy_taxfilt, file='Babiesphyfilt.RData')
