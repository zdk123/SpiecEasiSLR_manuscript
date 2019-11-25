## Make robPCA from Low Rank component ##
source('../../scripts/analyze_rob_pca_funs.R')
tmp <- build_pca('AGPNetFits.RData')
scores  <- tmp$scores
Rscores <- tmp$Rscores
dims    <- tmp$dims

## hardcode columns that should be factors (some look numeric)
col_fact <- read.table('map_fact_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]
## numeric columns
col_num <- read.table('map_num_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]

col_date <- c(collection_date='%m/%d/%Y')

map <- process_map(col_fact, col_num, col_date)
map$BarcodeSequence <- NULL

## remove NA cats
map <- map[,colMeans(apply(map, 2, is.na)) <= .4]
col_fact <- col_fact[col_fact %in% colnames(map)]
col_fact <- col_fact[apply(map[,col_fact], 2, function(x) length(unique(x)) < 50)]

map <- map[,c(col_fact[col_fact %in% colnames(map)],
              col_num[col_num %in% colnames(map)],
              col_date[col_date %in% colnames(map)])]

## attach depth data
map <- attach_depth(list.files(pattern='*filt.RData'), map)

save(scores, map, file='scores.RData')
## get all data
xCors <- cross_corr(scores, map, c(col_num, col_date), col_fact)

pvals <- signif_factors(scores, map)
pvals$n_features <- length(col_num) + length(col_fact) + length(col_date) + 2

keep  <- signif_features(pvals[[1]])
ord   <- compare_models(scores, Rscores, map[,keep[[1]]])
