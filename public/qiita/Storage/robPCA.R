## Make robPCA from Low Rank component ##
source('../../scripts/analyze_rob_pca_funs.R')
tmp <- build_pca('StorageNetFits.RData')
scores  <- tmp$scores
Rscores <- tmp$Rscores
dims    <- tmp$dims


## hardcode columns that should be factors (some look numeric)
col_fact <- read.table('map_fact_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]

map <- process_map(col_fact)
## attach depth data
map <- attach_depth(list.files(pattern='*filt.RData'), map)

## get all data
xCors <- cross_corr(scores, map, c(), col_fact)

pvals <- signif_factors(scores, map)
pvals$n_features <- length(col_num) + length(col_fact) + length(col_date) + 2

keep  <- signif_features(pvals[[1]])
ord   <- compare_models(scores, Rscores, map[,keep[[1]]])
