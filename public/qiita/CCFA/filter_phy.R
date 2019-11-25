library(phyloseq)
library(dplyr)
load('CCFAphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='CCFAphyfilt.RData')
