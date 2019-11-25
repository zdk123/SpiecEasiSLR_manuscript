library(phyloseq)
library(dplyr)
load('TwinsUKphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='TwinsUKphyfilt.RData')
