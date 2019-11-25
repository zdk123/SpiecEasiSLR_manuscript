library(phyloseq)

import_map <- function(x, ...) {
  tmp <- read.delim(x, comment.char="")
  rownames(tmp) <- tmp[,1]
  tmp
}


files  <- list.files('BIOM', pattern="*.biom$", recursive=TRUE, full.names=TRUE)
biomid <- basename(dirname(files))
mapfiles <- sapply(paste0("^", biomid , "_"), list.files, path='mapping_files', full.names=TRUE)

maps <- lapply(mapfiles, import_map)

phys <- lapply(files, import_biom)
for (i in 1:length(phys))
  sample_data(phys[[i]]) <- maps[[i]]

tax <- sapply(phys, taxa_names) %>% unlist %>% table
keeptax <- names(tax[tax>=length(phys)])

phy_filt <- lapply(phys, function(phy) prune_taxa(taxa_names(phy) %in% keeptax, phy))

phy <- do.call(merge_phyloseq, phy_filt)
save(phy, file='CCFAphy.RData')
