maxf1 <- new_metric("F1", "F1 score",
                metric = function(draw, model, out) {
                    max(out$pr$F1)
                })


hamfun <- function(x, y) {
  sum(x!=y)
}
hamming <- new_metric("hamming", "Hamming Distance",
              metric = function(draw, model, out) {
                min(out$vham)
              })

maxf1_2 <- new_metric("F1_2", "F1 score - 2",
                metric = function(draw, model, out) {
                    max(out$pr2$F1)
                })


hamming_2 <- new_metric("hamming_2", "Hamming Distance - 2",
              metric = function(draw, model, out) {
                min(out$vham2)
              })


zeronorm <- new_metric("Znorm", "Zero norm",
              metric = function(draw, model, out) {
                sum(out$true)/2
              })

zeronorm_2 <- new_metric("Znorm_2", "Zero norm -2",
              metric = function(draw, model, out) {
                sum(out$true2)/2
              })
