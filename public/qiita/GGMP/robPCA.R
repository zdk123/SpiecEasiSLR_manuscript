## Make robPCA from Low Rank component ##
source('../../scripts/analyze_rob_pca_funs.R')
tmp <- build_pca('GGMPNetFits.RData')
scores  <- tmp$scores
Rscores <- tmp$Rscores
dims    <- tmp$dims

## hardcode columns that should be factors (some look numeric)
col_fact <- read.table('map_fact_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]
## numeric columns
col_num <- read.table('map_num_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]


col_date <- c(dna_extraction_date='%Y.%m.%d')

map <- process_map(col_fact, col_num, col_date)
## attach depth data
map <- attach_depth(list.files(pattern='*filt.RData'), map)

## get all data
xCors <- cross_corr(scores, map, c(col_num, col_date), col_fact)


pvals <- signif_factors(scores, map)
pvals$n_features <- length(col_num) + length(col_fact) + length(col_date) + 2

keep  <- signif_features(pvals[[1]])
## remove NA cats
ord   <- compare_models(scores, Rscores, map[,keep[[1]]])


##
# glfun <- function(X, lambda, rho, p, ...) {
#   est <- QUIC::QUIC(cov(X, use='pairwise'), rho=rho, path=lambda, msg=0, tol=1e-2, ...)
#   path <- lapply(seq(length(lambda)), function(i) {
#               tmp <- est$X[,,i]; diag(tmp) <- 0
#               as(tmp!=0, "lgCMatrix")[1:p,1:p]
#   })
#   est$path <- path
#   est
# }
#
# ind <- c(86, 102, 147, 193, 229, 242, 272, 309, 313, 350, 382, 418, 472)
# subind <- which(ind %in% c(147, 382, 418))
# #load('ind.RData')
#
# Y <- model.matrix.lm(~., data=map[,keep[[1]]], na.action='na.pass')[,-1]
#
# #y.ind <- which(apply(abs(cor(tmp$X, Y, use='pairwise'))>=.2, 2, any))
#
# X.y <- cbind(tmp$X, scale(Y))
#
# #X.y <- t(na.omit(t(scale(na.omit(X.y)))))
# q <- ncol(X.y)
# p <- ncol(tmp$X)

#lammax <- max(abs(t(X.y) %*% X.y))/nrow(X.y)
#baserho <- rbind(matrix(0,q-p,q),cbind(matrix(1,p,p), matrix(0,p,q-p)))
# baserho <- c(rep(1,p), rep(0,q-p))
#
# args <- list(lambda=pulsar::getLamPath(1,1e-2,50,TRUE),
#              rho=baserho, p=p)
# fit <- pulsar::pulsar(X.y, glfun, args, thresh=.1, ncores=4)
# est <- pulsar::refit(fit)

#Theta <- fit$est$icov[[fit$stars$opt.index]]
#fit$stars$opt.index
#Theta <- est$est$X[,,50]
# Sig <- cov(X.y, use='pairwise.complete')
# Theta <- MASS::ginv(Sig)
#
# Theta_O  <- Theta[1:p,1:p]
# #Theta_O  <- Theta_O * sign(abs(S[ind,ind]))
# Theta_H  <- Theta[(p+1):q,(p+1):q]
# Theta_OH <- Theta[1:p,(p+1):q]
# K <- Theta_OH%*%MASS::ginv(Theta_H)%*%t(Theta_OH)

#Theta_OK <- solve(cor(tmp$X[,ind], use='pairwise.complete'))

#load('GGMPNetFits.RData')
# S <- se.slr[[4]]$est$icov[[se.slr[[4]]$select$stars$opt.index]]
# L <- se.slr[[4]]$est$resid[[se.slr[[4]]$select$stars$opt.index]]
# #Lsvd <- svd(L)
# #Y<- tmp$X %*% Lsvd$u[,1:9] %*% diag(1/Lsvd$d[1:9])
#
# # (Theta_O)[subind,subind]
# # (S)[ind,ind][subind,subind]
#
# metanet <- Theta[ind,(p+1):ncol(Theta)]
# metanet[abs(metanet)<=.5] <- 0
# metanet <- sign(abs(metanet)) * Sig[ind,(p+1):ncol(Theta)]
# metanet <- metanet[,y.ind2<-apply(abs(sign(metanet)), 2, sum)>0]
#
# ## combine cats
# colnames(metanet) <- stringr::str_match(colnames(metanet), "^([a-z_]*).*$")[,2]
# metanet <- sapply(unique(colnames(metanet)), function(j) rowSums(sign(abs(metanet[,colnames(metanet)==j,drop=FALSE]))))
# q <- p + ncol(metanet)
#
# augNet <- function(net, mnet=metanet) {
#   tmp <- rbind(cbind(net, mnet),
#                cbind(t(mnet), matrix(0, ncol(mnet), ncol(mnet))))
#   diag(tmp) <- 0
#   tmp
# }
#
# net <- as.matrix(cov2cor(Sig)[ind,ind] * sign(abs(S[ind,ind])))
# diag(net) <- 0
# net <- augNet(net)
#
# # net2 <- as.matrix(solve(S)[ind,ind]*sign(abs(S[ind,ind])))
# # diag(net2) <- 0
# # net2 <- augNet(net2)
#
# #net2 <- as.matrix(cov2cor(Sig)[ind,ind] * se.gl$refit$stars[ind,ind])
# S1 <- se.mb$refit$stars
# #net2 <- as.matrix(cov2cor(Sig)[ind,ind] * sign(abs(S1[ind,ind])))
# net2 <- as.matrix(cov2cor(MASS::ginv(K)))[ind,ind]
# diag(net2) <- 0
# #net2 <- net[1:length(ind),1:length(ind)]-net2
# #net2[abs(net2)<=.1] <- 0
#
# j <- c(rep(2,length(ind)), rep(1,q-p))
#
# par(mfrow=c(1,2))
# coord <- plot(network::network(net), usearrows=FALSE, edge.col=ifelse(net>0, 'forestgreen','red'), vertex.col=c(rep(2,length(ind)), rep(3,q-p)))
# plot(network::network(net2),  coord=coord, usearrows=FALSE, edge.col=ifelse(net2>0, 'forestgreen','red'))
#
# (Theta_O-K)[subind,subind]
# (S-L)[ind,ind][subind,subind]
