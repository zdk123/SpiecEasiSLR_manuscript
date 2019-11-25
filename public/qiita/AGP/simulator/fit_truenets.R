## define training models per method ##
true_ising <- function(X) {
  spiec.easi(X, method='ising', nlambda=200,
          lambda.min.ratio=1e-2, pulsar.params=list(ncores=30, thresh=0.05, rep.num=100),
          verbose=FALSE, lambda.log=FALSE)
}

true_poisson <- function(X) {
  spiec.easi(X, method='poisson', nlambda=200,
          lambda.min.ratio=1e-3, pulsar.params=list(ncores=30, thresh=0.05, rep.num=100),
          verbose=FALSE, lambda.log=FALSE)
}

true_glasso <- function(X) {
  spiec.easi(X, method='glasso', nlambda=200,
          lambda.min.ratio=1e-2, pulsar.params=list(ncores=30, thresh=0.05, rep.num=100),
          verbose=FALSE, lambda.log=FALSE)
}

true_mb <- function(X) {
  spiec.easi(X, method='mb', nlambda=200,
          lambda.min.ratio=1e-2, pulsar.params=list(ncores=30, thresh=0.05, rep.num=100),
          verbose=FALSE, lambda.log=FALSE)
}


true_slr <- function(X) {
  # r = 43 (double check)
  spiec.easi(X, method='slr', nlambda=200, r=r, lambda.max=1,
          lambda.min.ratio=1e-3, pulsar.params=list(ncores=30, thresh=0.05, rep.num=100),
          cor=TRUE, verbose=FALSE, lambda.log=FALSE)
}

common_scale <- function(x, depth) round(SpiecEasi::norm_to_total(x)*depth)

train_all_models <- function(X) {
  ## Binary data for ising
  Xs   <- sign(X)
  nobs <- colSums(Xs)
  n <- nrow(Xs)
  keep_s <- which(nobs>2 & nobs<(n-1))
  Xs <- Xs[,keep_s]

  ## common scale for poisson
  depths <- rowSums(X)
  Xc <- t(apply(X, 1, common_scale, depth=min(depths)))
  keep_c <- which(colSums(Xc)!=0)
  Xc <- Xc[,keep_c]

  print('mb')
  se.mb  <- true_mb(X)
  print('gl')
  se.gl  <- true_glasso(X)
  print('slr')
  se.slr <- true_slr(X)
  print('poi')
  se.poi <- true_poisson(Xc)
  print('is')
  se.is  <- true_ising(Xs)

  list(se.mb=(se.mb),
       se.gl=(se.gl),
       se.slr=(se.slr),
       se.poi=(se.poi),
       se.is=(se.is),
       keep_s=keep_s,
       keep_c=keep_c)
#       ind=ind)
}

true <- train_all_models(otu_train)

## alternative refits
minE <- min(sapply(true[1:5], function(x) sum(getRefit(x))))/2

filter_true <- function(true, nedges) {
  tmp <- getOptMerge(true)
  tmp[lower.tri(tmp)] <- 0
  set.seed(10010)
  tmp[upper.tri(tmp)] <- rank(SpiecEasi::triu(tmp, 1), ties.method='random')
  tmp <- forceSymmetric(tmp)
  tmp <- tmp <= minE
  diag(tmp) <- FALSE
  tmp
}
for (i in 1:5) {
  true[[i]]$refit2 <- filter_true(true[[i]], minE)
}

save(true, file='TrueNets.RData')
