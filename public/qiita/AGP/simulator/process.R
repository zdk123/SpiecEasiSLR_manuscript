
load('../AGPphyfilt.RData')

## define 'training' dataset as 9 largest runs
runtab <- sort(table(sample_data(phy_taxfilt)$run_prefix))
runs <- names(runtab)
train <- tail(runs, 8)
test <- setdiff(runs, train)


phy_train <- subset_samples(phy_taxfilt, run_prefix %in% train)
otu_train <- t(phy_train@otu_table@.Data)

phy_test <- subset_samples(phy_taxfilt, run_prefix %in% test)
otu_test <- t(phy_test@otu_table@.Data)

## Filter to small p, just for development purposes
# set.seed(10010)
# ind <- sample(taxa_names(phy_train), 30)
# otu_train <- otu_train[,ind]
# otu_test  <- otu_test[,ind]
