library(phyloseq)
library(dplyr)
load('Herfarthphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool")
save(phy_taxfilt, file='Herfarthphyfilt.RData')
