library(Matrix)


getOptSLR <- function(x) {
  ind <- x$ebic>1e-1
  x[ind][[which.min(x$ebic[ind])]]
}
getOpt <- function(x, sign=FALSE) {
  opt <- x
  if (is.null(opt$refit)) {
    opt <- getOptSLR(x)
  }
  gg_ids <- colnames(opt$est$data)
  if (is.null(gg_ids)) {
    gg_ids <- colnames(get('X', envir=opt$select$envir))
  }
  if ('stars' %in% names(opt$refit))
    tmp <- opt$refit$stars
  else
    tmp <- opt$refit

  colnames(tmp) <- rownames(tmp) <- gg_ids
  if (!sign) return(tmp)
  else {
    opt.index <- opt$select$stars$opt.index
    method <- eval(opt$select$call$fargs, opt$select$envir)$method
    method <- if(is.null(method)) 'slr' else method
    snet <- switch(
       method,
       glasso=,
       slr=sign(solve(opt$est$icov[[opt.index]])),
       poisson=, ising=,
       mb=sign(SpiecEasi:::symBeta(opt$est$beta[[opt.index]]))
      )
    return(tmp*snet)
  }
}

edge_density <- function(gr) {
    p <- ncol(gr)
    e <- sum(abs(gr))/2
    2*e/(p*(p-1))
}

proportion_pos <- function(gr) {
  pos <- sum(gr>0)
  pos / (sum(gr<0)+pos)
}

# mcc <- function(x, y) {
#   tmp <- common_tax(x, y)
#   x <- tmp$x
#   y <- tmp$y
#   p <- tmp$p
#
#   if (sum(x) < sum(y)) {
#     tmp <- x
#     x <- y
#     y <- tmp
#   }
#   p <- nrow(x)
#   maxe <- p*(p-1)/2
#   e <- sum(y)/2
#   tp <- sum((x==1)&(y==1))/(2*e)
#   tn <- sum((x==0)&(y==0))/sum(!y)
#   fp <- sum((x==1)&(y==0))/(2*(maxe-e))
#   fn <- sum((x==0)&(y==1))/(2*e)
#   num <- (tp*tn - fp*fn)
#   denom <- (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
#   num/sqrt(denom)
# }


f1 <- function(x, y) {
  #SpiecEasi::huge.pr(list(x), y, plot=FALSE, verbose=FALSE)$F1
  tmp <- common_tax(x, y)
  x <- tmp$x
  y <- tmp$y
  p <- tmp$p

  if (sum(x) < sum(y)) {
    tmp <- x
    x <- y
    y <- tmp
  }
  p <- nrow(x)
  maxe <- p*(p-1)/2
  e <- sum(y)/2
  tp <- sum((x==1)&(y==1))/(2*e)
  tn <- sum((x==0)&(y==0))/sum(!y)
  fp <- sum((x==1)&(y==0))/(2*(maxe-e))
  fn <- sum((x==0)&(y==1))/(2*e)
  2*tp / (2*tp + fp + fn)
}

# pairwise edge agreement
gr2edgelist <- function(obj) {
  adj <- as(obj, "symmetricMatrix")
  emat <- Matrix::summary(adj)[,1:2]
  cbind(rownames(adj)[emat[,1]], colnames(adj)[emat[,2]])
}

node_dist <- function(x, y) {
  x <- rownames(x)
  y <- rownames(y)
  1 - (length(intersect(x,y))/min(length(x), length(y)))
}

edge_dist2 <- function(x, y) {
  tmp <- common_tax(x, y)
  x <- gr2edgelist(abs(tmp$x))
  y <- gr2edgelist(abs(tmp$y))
  x <- apply(x, 1, paste, collapse=".")
  y <- apply(y, 1, paste, collapse=".")
  1 - (length(intersect(x,y))/min(length(x), length(y)))
}

edge_dist <- function(x,y) {
  x <- abs(x)
  y <- abs(y)
  ll <- sum((x==1)&(y==1))
  1 - (ll/min(sum(x), sum(y)))
}

## randomize edges in a graph
rgraph <- function(x) {
# x -> adjacency matrix
  suppressMessages(
  tmp <- SpiecEasi::triu2diag(sample(SpiecEasi::triu(x, k=1))))
  rownames(tmp) <- colnames(tmp) <- rownames(x)
  as(tmp, 'sparseMatrix')
}

common_tax <- function(x, y) {
  commontax <- intersect(rownames(x), rownames(y))
  p <- length(commontax)
  if (p > 0) {
    x <- x[commontax,commontax]
    y <- y[commontax,commontax]
  }
  list(x=x, y=y, p=p)
}

## function to extract/compare common subnetwork
X_dist <- function(x, y, f, ...) {
  tmp <- common_tax(x, y)
  x <- abs(tmp$x)
  y <- abs(tmp$y)
  p <- tmp$p
  dens <- 25
  if (p < 2*dens || (sum(x) < dens || sum(y) < dens)) return(NA)
  f(x, y, ...)
}

R_dist <- function(x, y, f, reps) {
  Reduce(function(x, y) mean(c(x,y), na.rm=TRUE),
    parallel::mclapply(1:reps, function(i) {
      X_dist(rgraph(x), rgraph(y), f=f)
   }, mc.cores=4)
  )
}


# netdist <- function(x, y, d='HIM') {
#   ## spectral distances ##
#   tmp <- common_tax(x, y)
#   x <- tmp$x
#   y <- tmp$y
#   suppressWarnings(
#   nettools::netdist(as.matrix(x), as.matrix(y),
#     d=d, components=FALSE, n.cores=1))
# }


## get net and attach OTU names ##
load_nets <- function(nets_rdata, ...) {
  print(nets_rdata)
  suppressWarnings(tmp <- load(nets_rdata))
  lapply(mget(tmp), getOpt, ...)
}

btw_cordist <- function(x, y) {
  tmp <- common_tax(x, y)
  x <- adj2igraph(abs(tmp$x))
  y <- adj2igraph(abs(tmp$y))
  cor(igraph::betweenness(x), igraph::betweenness(y), method='spearman')
}

degree_cordist <- function(x, y) {
  tmp <- common_tax(x, y)
  x <- abs(tmp$x)
  y <- abs(tmp$y)
  cor(rowSums(x), rowSums(y), method='spearman')
}

graphletcor_dist <- function(x, y) {
  # tmp <- common_tax(x, y)
  x <- abs(x)
  y <- abs(y)
  dist(rbind(pulsar::gcvec(x), pulsar::gcvec(y)))
}
