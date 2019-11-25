library(phyloseq)
library(dplyr)
load('Familyphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='Familyphyfilt.RData')
