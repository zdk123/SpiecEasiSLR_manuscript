library(phyloseq)
library(dplyr)
source('scripts/filter_fun.R')
currdir <- getwd()

projects <- sapply(list.dirs(recursive=FALSE),
              list.dirs, recursive=FALSE) %>%
              unlist %>% unname
dimmat <- matrix(NA, nrow=length(projects), ncol=4)
rownames(dimmat) <- projects
colnames(dimmat) <- c('P', 'N', 'p', 'n')

for (d in projects) {
  setwd(d)
  print(d)
  source('filter_phy.R')
  print(phy_taxfilt)
  dimmat[d,1] <- ntaxa(phy)
  dimmat[d,2] <- nsamples(phy)
  dimmat[d,3] <- ntaxa(phy_taxfilt)
  dimmat[d,4] <- nsamples(phy_taxfilt)
  rm(phy, phy_taxfilt)
  setwd(currdir)
}
