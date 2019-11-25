## @knitr metrics

hamming <- new_metric("hamming", "hamming distance",
                 metric = function(draw, model, out) {
                    out$minHam
                 })

maxf1 <- new_metric("F1", "F1 score",
                 metric = function(draw, model, out) {
                    out$maxF1
                 })

perfectf1 <- new_metric("perfectF1", "perfect F1 score",
                 metric = function(draw, model, out) {
                    out$maxF1 == 1
                 })

oracle_size <- new_metric("oracle_znorm", "oracle znorm",
                 metric = function(draw, model, out) {
                   sum(out$sel2)/2
                 })

zeronorm.cov <-  new_metric("Zero_cov", "Zero norm cov",
                 metric = function(draw, model, out) {
                    .5*(sum(abs(sign(draw$Cov))) - ncol(draw$Cov))
                 })

zeronorm.icov <-  new_metric("Zero_icov", "Zero norm icov",
                 metric = function(draw, model, out) {
                   graph <- abs(sign(out$true))
                   diag(graph) <- 0
                    sum(graph)/2
                 })

varrowmeans.icov <-  new_metric("varrowmeans_icov", "Variance of Rowmeans icov",
                metric = function(draw, model, out) {
                   var(rowMeans(draw$Prec))
                })


maxdeg <-  new_metric("maxdeg", "Maximum degree",
                 metric = function(draw, model, out) {
                    max(rowSums(draw$graph))
                 })

# d <-  new_metric("d", "degree",
#                  metric = function(draw, model, out) {
#                     max(rowSums(draw$graph))
#                  })

spread.incoh <-  new_metric("spread_incoh", "spread / incoherence",
                 metric = function(draw, model, out) {
                    draw$spread / draw$incoh
                 })


incoh <-  new_metric("incoh", "incoherence",
                 metric = function(draw, model, out) {
                    draw$incoh
                 })

incoh2 <-  new_metric("incoh2", "incoherence2",
                metric = function(draw, model, out) {
                   draw$incoh2
                })

mu <-  new_metric("mu", "mu",
                 metric = function(draw, model, out) {
                    draw$mu
                 })

eta <-  new_metric("eta", "eta",
                 metric = function(draw, model, out) {
                    draw$eta
                 })


spread <-  new_metric("spread", "spread",
                 metric = function(draw, model, out) {
                    draw$spread
                 })
