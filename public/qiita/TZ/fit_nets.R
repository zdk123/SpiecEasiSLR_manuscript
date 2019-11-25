library(Matrix)
library(phyloseq)
library(SpiecEasi)

load('TZphyfilt.RData')

#load('TZNetFits.RData')


se.gl <- spiec.easi(phy_taxfilt, method='glasso', nlambda=50,
            lambda.min.ratio=1e-1, lambda.log=FALSE,
            pulsar.params=list(ncores=38, rep.num=35))

se.mb <- spiec.easi(phy_taxfilt, method='mb', nlambda=50,
            lambda.min.ratio=1e-1, lambda.log=FALSE,
            pulsar.params=list(ncores=38, rep.num=35))

X <- sign(phy_taxfilt@otu_table@.Data)
phy_tmp <- prune_taxa(apply(X, 1, var) != 0, phy_taxfilt)
se.is <- spiec.easi(phy_tmp, method='ising', nlambda=50,
            lambda.min.ratio=1e-1, lambda.log=FALSE,
            pulsar.params=list(ncores=38, rep.num=35))

se.poi <- spiec.easi(phy_taxfilt, method='poisson', nlambda=50,
            lambda.min.ratio=1e-2, lambda.log=FALSE,
            pulsar.params=list(ncores=35, rep.num=35))

#save(se.gl, se.mb, se.is, se.poi, file='NetFits.RData')
ranks <- round(exp(seq(log(2), log(32), len=6)))
se.slr <- spiec.easi(phy_taxfilt, method='slr', nlambda=50,
            lambda.min.ratio=1e-2, r=ranks, lambda.log=TRUE,
            pulsar.params=list(ncores=32, rep.num=30))
se.slr$ebic <- sapply(se.slr, function(x)
                        ebic(x$refit$stars, x$est$data,
                             x$est$loglik[x$select$stars$opt.index]))

save(se.gl, se.mb, se.is, se.poi, se.slr, file='TZNetFits.RData')
