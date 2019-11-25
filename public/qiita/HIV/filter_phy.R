library(phyloseq)
library(dplyr)
load('HIVphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='HIVphyfilt.RData')
