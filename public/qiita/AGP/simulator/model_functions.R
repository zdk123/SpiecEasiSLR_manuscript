## call from main.R
make_AGP_sub <- function(n) {
  new_model(name = "AGP",
            label = sprintf("n = %s", n),
            params = list(n=n),
            simulate = function(n, nsim) {
              if (n>nrow(otu_train)) stop('max sample size exceeded')
              pcoef <- n/sum(runtab[test])

              lapply(1:nsim, function(i) {
                  ## subset equal proportions of each dataset
                  subsamp <- sample(1:nsamples(phy_test), n)
                  list(X=otu_test[subsamp,])
                })
            }
  )
}
