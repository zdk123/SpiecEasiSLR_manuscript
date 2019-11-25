library(phyloseq)
library(dplyr)
load('Mayophy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='Mayophyfilt.RData')
