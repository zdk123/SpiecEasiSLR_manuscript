hamfun <- function(x, y) {
  sum(x!=y)
}

se.slr <- new_method("se-slr", "SPIEC-EASI SLR",
        method = function(model, draw) {
          X <- draw$X+1
          est <- sparseLowRankiCov(cor(t(clr(X,1))), lambda.min.ratio=1e-3, lambda.max=1.5, nlambda=100, r=r, cor=TRUE)
          est$keep <- 1:ncol(X)
          est$fun <- sparseLowRankiCov
          est$true <- getRefit(true$se.slr)
          est$true2 <- true$se.slr$refit2
          est$optmerge <- getOptMerge(true$se.slr)
          est$pr   <- huge.pr(est$path, est$true, FALSE, FALSE)
          est$vham <- sapply(est$path, hamfun, y=est$true)
          est$pr2   <- huge.pr(est$path, est$true2, FALSE, FALSE)
          est$vham2 <- sapply(est$path, hamfun, y=est$true2)
          return(est)
        }
      )


se.gl <- new_method("se-gl", "SPIEC-EASI GLASSO",
        method = function(model, draw) {
          X <- draw$X+1
          est <- sparseiCov(cor(t(clr(X,1))), method='glasso', lambda.min.ratio=5e-3, nlambda=100)
          est$keep <- 1:ncol(X)
          est$fun <- sparseiCov
          est$true <- getRefit(true$se.gl)
          est$true2 <- true$se.gl$refit2
          est$optmerge <- getOptMerge(true$se.gl)
          est$pr <- huge.pr(est$path, est$true, FALSE, FALSE)
          est$vham <- sapply(est$path, hamfun, y=est$true)
          est$pr2   <- huge.pr(est$path, est$true2, FALSE, FALSE)
          est$vham2 <- sapply(est$path, hamfun, y=est$true2)
          return(est)
        }
      )

se.mb <- new_method("se-mb", "SPIEC-EASI MB",
        method = function(model, draw) {
          X <- draw$X+1
          est <- sparseiCov(cor(t(clr(X,1))), method='mb',
                            lambda.min.ratio=5e-3, nlambda=100)
          est$keep <- 1:ncol(X)
          est$fun <- sparseiCov
          est$true <- getRefit(true$se.mb)
          est$true2 <- true$se.mb$refit2
          est$optmerge <- getOptMerge(true$se.mb)
          est$pr <- huge.pr(est$path, est$true, FALSE, FALSE)
          est$vham <- sapply(est$path, hamfun, y=est$true)
          est$pr2   <- huge.pr(est$path, est$true2, FALSE, FALSE)
          est$vham2 <- sapply(est$path, hamfun, y=est$true2)
          return(est)
        }
      )

se.is <- new_method("se-is", "SPIEC-EASI ising",
        method = function(model, draw) {
          X <- draw$X
          Xs   <- sign(X)[,true$keep_s]
          nobs <- colSums(Xs)
          n <- nrow(Xs)
          keep_s <- which(nobs>2 & nobs<(n-1))
          est <- spiec.easi(Xs[,keep_s], method='ising',
                      nlambda=100, lambda.min.ratio=1e-3,
                      verbose=FALSE, pulsar.select=FALSE)$est
          est$keep <- keep_s
          est$true <- getRefit(true$se.is)[keep_s,keep_s]
          est$true2 <- true$se.is$refit2[keep_s,keep_s]
          est$optmerge <- getOptMerge(true$se.is)
          est$pr <- huge.pr(est$path, est$true, FALSE, FALSE)
          est$vham <- sapply(est$path, hamfun, y=est$true)
          est$pr2   <- huge.pr(est$path, est$true2, FALSE, FALSE)
          est$vham2 <- sapply(est$path, hamfun, y=est$true2)
          return(est)
        }
      )

common_scale <- function(x, depth) round(SpiecEasi::norm_to_total(x)*depth)
se.poi <- new_method("se-poi", "SPIEC-EASI poisson",
        method = function(model, draw) {
          X <- draw$X
          depths <- rowSums(X)
          Xc <- t(apply(X, 1, common_scale, depth=min(depths)))[,true$keep_c]
          keep_c <- which(colSums(Xc)!=0)
          est <- spiec.easi(Xc[,keep_c], method='poisson', nlambda=100,
                  lambda.min.ratio=1e-3, verbose=FALSE,
                  pulsar.select=FALSE)$est
          est$keep <- keep_c
          est$true  <- getRefit(true$se.poi)[keep_c,keep_c]
          est$true2 <- true$se.poi$refit2[keep_c,keep_c]
          est$optmerge <- getOptMerge(true$se.poi)
          est$pr <- huge.pr(est$path, est$true, FALSE, FALSE)
          est$vham <- sapply(est$path, hamfun, y=est$true)
          est$pr2   <- huge.pr(est$path, est$true2, FALSE, FALSE)
          est$vham2 <- sapply(est$path, hamfun, y=est$true2)
          return(est)
        }
      )
