library(dplyr)
library(xtable)
meta <- read.table('data/metadata.txt', sep="\t", header=TRUE)

meta %>%
  select(name, qiita_ids, gene, region, sequencing_platform,
         gut_samples_filt, OTUs, OTUs_filt, citekey) %>%
  mutate(citekey=paste0("\\cite{", citekey, "}")) -> meta
