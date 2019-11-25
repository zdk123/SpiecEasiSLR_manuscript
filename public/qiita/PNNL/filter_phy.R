library(phyloseq)
library(dplyr)
load('PNNLphy.RData')

phy_taxfilt <- filter(phy, 'sample_type', "stool",
#  getmindepth=function(depth) quantile(depth, .7),
  gettaxind=function(phy) rowSums(sign(phy@otu_table@.Data)) > nsamples(phy)*.6)


save(phy_taxfilt, file='PNNLphyfilt.RData')
