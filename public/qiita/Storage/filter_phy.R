library(phyloseq)
library(dplyr)
load('Storagephy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='Storagephyfilt.RData')
