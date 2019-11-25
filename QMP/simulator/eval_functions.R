corshrink <- function(S) {
  as.matrix(corpcor::cor.shrink(cov2cor(S), verbose=FALSE))
}

covshrink <- function(S) {
  as.matrix(corpcor::cov.shrink((S), verbose=FALSE))
}

f <-  cov2cor #force #
fshrink <- corshrink # covshrink #

minFNorm <- new_metric("Fnorm", "Fnorm",
                metric = function(draw, model, out) {
                  S <- f(draw$Sig)
                  Fnorm <- sapply(out$cov, function(x) norm(as.matrix(f(x))-S, 'F'))
                  min(Fnorm)
                 })

oracleNorm <- new_metric("OracleFnorm", "Oracle Fnorm",
                metric = function(draw, model, out) {
                  S <- f(draw$Sig)
                  Fnorm <- sapply(out$cov, function(x) norm(as.matrix(f(x))-S, 'F'))
                  norm(f(out$cov[[which.min(Fnorm)]]), 'F')
                })


SigNorm <- new_metric("SigNorm", "Sigma Fnorm",
                metric = function(draw, model, out) {
                  S <- f(draw$Sig)
                  norm(S, 'F')
                })


minFNorm2 <- new_metric("Fnorm_shrink", "Fnorm shrink",
                metric = function(draw, model, out) {
                  S <- fshrink(f(draw$Sig))
                  Fnorm <- sapply(out$cov, function(x) norm(as.matrix(f(x))-S, 'F'))
                  min(Fnorm)
                })

oracleNorm2 <- new_metric("OracleFnorm_shrink", "Oracle Fnorm Shrink",
                metric = function(draw, model, out) {
                  S <- fshrink(f(draw$Sig))
                  Fnorm <- sapply(out$cov, function(x) norm(as.matrix(f(x))-S, 'F'))
                  norm(f(out$cov[[which.min(Fnorm)]]), 'F')
                })

SigNorm2 <- new_metric("SigNorm_shrink", "Sigma Fnorm Shrink",
                metric = function(draw, model, out) {
                  S <- fshrink(f(draw$Sig))
                  norm(S, 'F')
                })


minFNorm3 <- new_metric("Fnorm_fshrink", "Fnorm Cov shrink",
                metric = function(draw, model, out) {
                  S <- fshrink(draw$Sig)
                  Fnorm <- sapply(out$cov, function(x) norm(as.matrix(cov2cor(x))-S, 'F'))
                  min(Fnorm)
                })

oracleNorm3 <- new_metric("OracleFnorm_fshrink", "Oracle Fnorm Cov Shrink",
                metric = function(draw, model, out) {
                  S <- fshrink(draw$Sig)
                  Fnorm <- sapply(out$cov, function(x) norm(as.matrix(cov2cor(x))-S, 'F'))
                  norm(f(out$cov[[which.min(Fnorm)]]), 'F')
                })

SigNorm3 <- new_metric("SigNorm_fshrink", "Sigma Fnorm Cov Shrink",
                metric = function(draw, model, out) {
                  S <- fshrink(draw$Sig)
                  norm(S, 'F')
                })

# incoh <-  new_metric("incoh", "incoherence",
#                 metric = function(draw, model, out) {
#                    draw$incoh
#                 })
#
# incoh2 <-  new_metric("incoh2", "incoherence2",
#                metric = function(draw, model, out) {
#                   draw$incoh2
#                })
#
# mu <-  new_metric("mu", "mu",
#                 metric = function(draw, model, out) {
#                    draw$mu
#                 })
#
# eta <-  new_metric("eta", "eta",
#                 metric = function(draw, model, out) {
#                    draw$eta
#                 })
#
#
# spread <-  new_metric("spread", "spread",
#                 metric = function(draw, model, out) {
#                    draw$spread
#                 })
