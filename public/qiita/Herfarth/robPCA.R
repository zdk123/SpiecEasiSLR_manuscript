## Make robPCA from Low Rank component ##
source('../../scripts/analyze_rob_pca_funs.R')
tmp <- build_pca('HerfarthNetFits.RData')
scores  <- tmp$scores
Rscores <- tmp$Rscores
dims    <- tmp$dims

## hardcode columns that should be factors (some look numeric)
col_fact <- read.table('map_fact_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]
## numeric columns
col_num <- read.table('map_num_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]


col_date <- c(collection_date="%m/%d/%y", collection_time="%H:%M", collection_timestamp="%m/%d/%y %H:%M")

map <- process_map(col_fact, col_num, col_date)
map$extraction_robot <- NULL
map$BarcodeSequence <- NULL
## attach depth data
map <- attach_depth(list.files(pattern='*filt.RData'), map)

## get all data
xCors <- cross_corr(scores, map, c(col_num, col_date), col_fact)

pvals <- signif_factors(scores, map)
pvals$n_features <- length(col_num) + length(col_fact) + length(col_date) + 2
keep  <- signif_features(pvals[[1]])
## remove NA cats
ord   <- compare_models(scores, Rscores, map[,keep[[1]]])
