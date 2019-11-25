library(phyloseq)
library(dplyr)

currdir <- getwd()

projects <- sapply(list.dirs(recursive=FALSE),
              list.dirs, recursive=FALSE) %>%
              unlist %>% unname

for (d in projects {
  setwd(d)
  print(d)
  source('fit_nets.R')
  rm(se.gl, se.mb, se.is, se.poi, se.slr)
  gc()
  setwd(currdir)
}
