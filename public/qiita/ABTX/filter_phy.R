library(phyloseq)
library(dplyr)
load('ABTXphy.RData')
# source('../../filter_fun.R')

phy_taxfilt <- filter(phy, 'env_material', "feces",
#  getmindepth=function(depth) quantile(depth, .7),
  gettaxind=function(phy) rowSums(sign(phy@otu_table@.Data)) > nsamples(phy)*.6)


save(phy_taxfilt, file='ABTXphyfilt.RData')
