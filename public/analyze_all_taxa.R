library(phyloseq)

taxa <- lapply(list.files(pattern="*phy.RData", recursive=TRUE), function(pdata) {
  load(pdata)
  taxa_names(phy)
})

taxafilt <- lapply(list.files(pattern="*phyfilt.RData", recursive=TRUE), function(pdata) {
  load(pdata)
  taxa_names(phy_taxfilt)
})

length(unique(unlist(taxa)))
length(unique(unlist(taxafilt)))


nsamples <- sapply(list.files(pattern="*phy.RData", recursive=TRUE), function(pdata) {
  load(pdata)
  nsamples(phy)
})

nsamples_filt <- sapply(list.files(pattern="*phyfilt.RData", recursive=TRUE), function(pdata) {
  load(pdata)
  nsamples(phy_taxfilt)
})
