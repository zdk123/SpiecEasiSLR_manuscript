library(phyloseq)
library(dplyr)
load('Bangladeshphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='Bangladeshphyfilt.RData')
