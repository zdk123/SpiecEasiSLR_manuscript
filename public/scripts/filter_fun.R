library(dplyr)

filttaxind   <- function(phy)
  rowSums(sign(phy@otu_table@.Data)) > nsamples(phy)*.37

filtmindepth <- function(depth) quantile(depth, .1)


filter <- function(phy, col='body_site', pattern="UBERON:feces",
                   getmindepth=filtmindepth, gettaxind=filttaxind) {
  ## Filter gut samples ##
  keep <- grepl(pattern, sample_data(phy)[,col]@.Data[[1]])
  phy_filt <- prune_samples(keep, phy)

  ## Filter low abundance samples ##
  depth <- colSums(phy_filt@otu_table@.Data)
  phy_filt2 <- prune_samples(depth > getmindepth(depth), phy_filt)

  phy_taxfilt <- phy_filt2 %>% prune_taxa(gettaxind(.), .)
  invisible(phy_taxfilt)
}
