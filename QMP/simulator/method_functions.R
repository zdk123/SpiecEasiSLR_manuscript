eps <- 1/10 #
#dname <- 'otus'
# eps <- 0
dname <- 'qmp'

se.slr <- new_method("seslr", "SPIEC-EASI SLR",
        method = function(model, draw) {
          X <- draw[[dname]]+eps
          est <- sparseLowRankiCov(cor(t(clr(X,1))), lambda.min.ratio=1e-2, lambda.max=2, nlambda=50, r=2)
          est$fun <- sparseLowRankiCov
          est$cov <- lapply(est$icov, function(x) forceSymmetric(solve(x)))
          return(est)
        }
    )

se.gl <- new_method("segl", "SPIEC-EASI glasso",
        method = function(model, draw) {
          X <- draw[[dname]]+eps
          est <- sparseiCov(cor(t(clr(X,1))), method='glasso', lambda.min.ratio=1e-2, nlambda=50)
          est$fun <- sparseiCov
          return(est)
        }
    )

coat <- new_method("coat", "COAT",
        method = function(model, draw) {
          X <- draw[[dname]]+eps
          est <- SpiecEasi::coat(cor(t(clr(X,1))), lambda.min.ratio=1e-2, nlambda=50, adaptive=FALSE)
          est$fun <- coat
          return(est)
        }
    )

crude <- new_method("crude", "CRUDE",
        method = function(model, draw) {
          X <- draw[[dname]]+eps
          est <- SpiecEasi::coat(cor(t(apply(X,1,SpiecEasi::norm_to_total))), lambda.min.ratio=1e-2, nlambda=50, adaptive=FALSE)
          est$fun <- coat
          return(est)
        }
    )
