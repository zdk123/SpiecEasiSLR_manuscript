library(phyloseq)
library(dplyr)
source('scripts/process_fun.R')

currdir <- getwd()

projects <- sapply(list.dirs(recursive=FALSE),
              list.dirs, recursive=FALSE) %>%
            unlist %>% unname

dimmat <- matrix(NA, nrow=length(projects), ncol=2)
rownames(dimmat) <- projects
taxa <- vector('list', length(projects))
names(taxa) <- projects

for (d in projects) {
  setwd(d)
  print(d)
  phy <- process(d)
  dimmat[d,1] <- ntaxa(phy)
  dimmat[d,2] <- nsamples(phy)
  taxa[[d]] <- taxa_names(phy)
  rm(phy)
  setwd(currdir)
}
