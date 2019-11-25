## @knitr methods
coat <- new_method("coat", "COAT",
          method = function(model, draw) {
            S <- draw$S
            est <- SpiecEasi::coat(S, nlambda=100, lambda.min.ratio=1e-3, adaptive=FALSE, shrinkDiag=FALSE, ret.icov=FALSE)
#            est$path <- lapply(est$cov, function(cov) as(abs(MASS::ginv(as.matrix(cov))) > 1e-4, 'lsCMatrix'))
#            draw$graph <- abs(cov2cor(draw$Cov))>.5
#            diag(draw$graph) <- 0
#            est$path <- lapply(est$icov, function(x) {tmp <- abs(sign(x)) ; diag(tmp) <- 0 ; tmp})
            pr  <- huge.pr(est$path, draw$graph, verbose=FALSE, plot=FALSE)
            ind <- which.max(pr$F1)
            est$sel.omega <- est$icov[[ind]]
            est$sel <- est$path[[ind]]

            ham <- sapply(est$path, function(x) sum(x!=as.matrix(draw$graph))/2)
            ind2 <- which.min(ham)
            est$sel2 <- est$path[[ind2]]
            est$minHam <- ham[ind2]

            est$true <- draw$Cov*draw$graph #draw$Prec #
            diag(est$true) <- diag(draw$Cov)
            est$maxF1 <- pr$F1[ind]
            est$pr <- pr
            return(est)
          }
    )

se.slr <- new_method("seslr", "SPIEC-EASI SLR",
          method = function(model, draw) {
            S <- draw$S
#            max  <- pulsar::getMaxCov(S, diag=TRUE)
#            lams <- pulsar::getLamPath(max, max*1e-2, 10)
            est <-  sparseLowRankiCov(S, r=2, nlambda=100, lambda.min.ratio=1e-3, lambda.max=.8, cor=FALSE)
            pr  <- huge.pr(est$path, draw$graph, verbose=FALSE, plot=FALSE)
            ind <- which.max(pr$F1)
            est$sel <- est$path[[ind]]
            est$sel.omega <- est$icov[[ind]]

            ham <- sapply(est$path, function(x) sum(x!=as.matrix(draw$graph))/2)
            ind2 <- which.min(ham)
            est$sel2 <- est$path[[ind2]]
            est$minHam <- ham[ind2]

            est$true <- draw$Prec
            est$maxF1 <- pr$F1[ind]
            est$pr <- pr
            return(est)
          }
    )

se.gl <- new_method("segl", "SPIEC-EASI glasso",
          method = function(model, draw) {
            S <- draw$S
#            max  <- pulsar::getMaxCov(S, diag=TRUE)
#            lams <- pulsar::getLamPath(max, max*1e-2, 10)
            est <- sparseiCov(S, method='glasso', lambda.min.ratio=1e-3, nlambda=100)
            pr  <- huge.pr(est$path, draw$graph, verbose=FALSE, plot=FALSE)
            ind <- which.max(pr$F1)
            est$sel <- est$path[[ind]]
            est$sel.omega <- est$icov[[ind]]

            ham <- sapply(est$path, function(x) sum(x!=as.matrix(draw$graph))/2)
            ind2 <- which.min(ham)
            est$sel2 <- est$path[[ind2]]
            est$minHam <- ham[ind2]

            est$true <- draw$Prec
            est$maxF1 <- pr$F1[ind]
            est$pr <- pr
            return(est)
          }
    )
