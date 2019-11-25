#!/bin/env Rscript
## random graphs p=16, p=102
## 3 or 4 latent effects
## Show SE glasso vs slr results:


library(SpiecEasi)

ind <- c(1:5, 13, 7:8, 11:12, 14, 6, 9:10, 15:16)

set.seed(10010)
g1 <- make_graph("erdos_renyi", 16, 10)
p1 <- graph2prec(g1, targetCondition=10, posThetaLims=c(1,1.5))
diag(p1) <- diag(p1) + rnorm(16, 0, .4)
#c1 <- prec2cov(p1)
## Add 'diet' block
p2 <- c(rep(.65, 5), rep(0, 16-5))
#p2 <- c(rep(0, 8), rep(.65, 2), rep(0, 4), rep(.65, 2))
p2 <- p2 %*% t(p2)
## add compositional effects
G <- diag(16)-1/16
p12 <- p1-p2
p12 <- Matrix::nearPD(p12)$mat
p3 <- G%*%(p12)%*%G - p12

j <- rep(1, 16)/16
p4 <- p12+p3
p4 <- p4[ind, ind]
g1 <- g1[ind,ind]

CCov <- as.matrix(G%*%solve(p12)%*%G)[ind, ind]
ACov <- as.matrix(    solve(p12)    )[ind, ind]
slr1 <- sparseLowRankiCov(as.matrix((CCov)), r=1, nlambda=100, lambda.min.ratio=1e-2)

slr <- sparseLowRankiCov(as.matrix((CCov)), r=2, nlambda=100, lambda.min.ratio=1e-2)
#gl <- sparseLowRankiCov(as.matrix((CCov)), r=0, nlambda=100, lambda.min.ratio=1e-5)
gl  <- sparseiCov(as.matrix((CCov)), method='glasso', lambda.min.ratio=1e-2, nlambda=150)

print(slr1.roc <- huge::huge.roc(slr1$path, as.matrix(g1)))
print(slr.roc <- huge::huge.roc(slr$path, as.matrix(g1)))
print(gl.roc  <- huge::huge.roc(gl$path, as.matrix(g1)))
dev.off()
slr.opt <- slr$icov[[which.max(slr.roc$F1)]]
resid  <- slr$resid[[which.max(slr.roc$F1)]]
gl.opt <- gl$icov[[which.max(gl.roc$F1)]]

library(ggplot2)
range0 <- function(x) {
  ## does range of x overlap 0
  xrange <- range(x)
  min(xrange)<0 & max(xrange)>0
}

## use a common scale for coloring purposes:
slr2 <- as.matrix(slr.opt-resid)

optrange <- c(min(min(gl.opt), min(slr.opt), min(slr2)),
              max(max(gl.opt), max(slr.opt), max(slr2)))
slr.opt <- as.matrix(scales::rescale(as.matrix(slr.opt), optrange, range(slr.opt))*abs(sign(slr.opt)))
gl.opt <- as.matrix(scales::rescale(as.matrix(gl.opt), optrange, range(gl.opt))*abs(sign(gl.opt)))
slr2 <- as.matrix(scales::rescale(as.matrix(slr2), optrange, range(slr2))*abs(sign(slr2)))


ggImage <- function(x, diag=TRUE, min1=min(slr2),
                    max1=max(slr2-diag(diag(slr2))), max2=max(slr2)) {
  x <- as.matrix(x)
  xdf <- reshape2::melt(x)
  p1 <- ggplot(aes(Var1, Var2), data=xdf) +
    geom_tile(aes(fill=value), color='grey70') +
    theme_classic() +
    xlab("") + ylab("") +
    coord_fixed(expand=FALSE) +
    scale_y_reverse(breaks=seq.int(ncol(x))) +
    scale_x_continuous(breaks=seq.int(nrow(x))) +
    theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

    if (diag) {
        cols <- c("red", 'white', "dodgerblue", 'dodgerblue4')
        vals <- c(min1,
                  ifelse(range0(x), 0, mean(x)),
                  max1, max2)
        rangefrom <- c(min1, max2)
    } else {
        if (all(range(x) < 0)) {
          cols <- c("red", 'white')
          vals <- c(min1, 0)
          rangefrom <- c(min(x), 0)
        } else if (all(range(x) > 0)) {
          cols <- c('white', 'dodgerblue4')
          vals <- c(0, max1)
          rangefrom <- c(0, max(x))
        } else {
          cols <- c("red", 'white', "dodgerblue")
          vals <- c(min1, ifelse(range0(x), 0, mean(x)), max1)
          rangefrom <- range(x) #c(min1, max1)
        }
    }
    vals <- scales::rescale(vals, from=rangefrom)
    p1 + scale_fill_gradientn(colors=cols, values=vals)
}

ggImage2 <- function(x) {
  x <- as.matrix(x)
  xdf <- reshape2::melt(x)
  ggplot(aes(Var1, Var2), data=xdf) +
    geom_tile(aes(fill=value), color='grey70') +
    theme_classic() +
    xlab("") + ylab("") +
    coord_fixed(expand=FALSE) +
    scale_y_reverse(breaks=seq.int(ncol(x))) +
    scale_x_continuous(breaks=seq.int(nrow(x))) +
    theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
}

rsvd <- svd(resid)
r2 <- (-rsvd$v[,2] %*% t(rsvd$v[,2]))
r1 <- (-rsvd$v[,1] %*% t(rsvd$v[,1]))

#r2 <- scales::rescale(r2, range(resid), range(r2))

set.seed(10010)
X <- exp(rmvnorm(20, rep(0,16), cov2cor(solve(p12))))
Z <- t(clr(X, 1))
loadings <- sqrt(diag(1/rsvd$d[1:2])) %*% t(rsvd$v[,1:2])
rob_scores <- scale(scale(Z) %*% t(loadings))
## fake data classes
cl <- (scale(Z)[,1] > -.16 & Z[,3] > -.3) + 1
library(xkcd)
scoresdf <- data.frame(rob_scores, class=as.factor(cl))


pdf('compositional.pdf', width=11, height=8)
gridExtra::grid.arrange(ggImage(slr2),
                        ggImage(gl.opt),
                        ggImage(slr.opt),
                        ggImage(r2, FALSE),
                        ggImage(r1, FALSE),
#                        ggImage(-resid),
                        ggImage( rsvd$v[,2,drop=FALSE]/2, FALSE),
                        ggImage(-rsvd$u[,2,drop=FALSE]/2, FALSE),
                        ggImage(rsvd$v[,c(1),drop=FALSE]/2, FALSE),
                        ggImage(loadings, FALSE),
                        ncol=3)

gridExtra::grid.arrange(
    ggImage2(t(Z[order(cl),])) + scale_fill_distiller(type='div', palette=4),
    ggImage2(rob_scores[order(cl),]) +
        scale_fill_distiller(type='div', palette=4),
    ncol=2
  )

dev.off()

pdf("compositional2.pdf", height=4, width=4)
plot(ggplot(aes(X1, X2), data=scoresdf) +
    geom_point(aes(color=class), size=4) +
    xkcdaxis(ceiling(range(rob_scores[,1])*5)/5, ceiling(range(rob_scores[,2])*10)/10) +
    guides(color=FALSE) +
    scale_color_brewer(palette='Set1') +
    xlab("PC1") + ylab("PC2"))
dev.off()
