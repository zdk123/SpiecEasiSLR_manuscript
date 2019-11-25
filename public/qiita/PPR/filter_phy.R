library(phyloseq)
library(dplyr)
load('PPRphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='PPRphyfilt.RData')
