ppareto <- function(x, scale, shape)
{
  ifelse(x > scale, 1 - (scale / x) ^ shape, 0)
}

qpareto <- function(y, scale, shape)
{
  ifelse(
    y >= 0 & y <= 1,
    scale * ((1 - y) ^ (-1 / shape)),
    NaN
  )
}


rand_vect <- function(N, M, maxdeg) {
#  scale <- 1
#  shape <- .5
#  lower_bound <- 1
#  upper_bound <- maxdeg+1

#  (quantiles <- ppareto(c(lower_bound, upper_bound), scale, shape))
#  uniform_random_numbers <- runif(N, quantiles[1], quantiles[2])
#  tmp <- round(qpareto(uniform_random_numbers, scale, shape))-1
  ind <- 0:maxdeg
  prob <-norm_to_total(ppareto(rev(ind+2), 1, 1))
  names(prob) <- ind
  tmp <- sample(ind, N, replace=TRUE, prob=prob)
  if ((sum(tmp) %% 2) != 0) {
    ind <- which.min(tmp)
    tmp[ind] <- tmp[ind]+1
  }
  tmp
}

connected.components <- function(A) {
  memb <- igraph::components(SpiecEasi::adj2igraph(A),
                            'strong')$membership
  tmp <- sapply(unique(memb), function(i) memb==i)
  tmp[,colSums(tmp)>1,drop=FALSE]
}

bridge_edges <- function(graph, exclude_degree_1=FALSE) {
  igr <- adj2igraph(graph)
  bicc <- igraph::biconnected.components(igr)
  ## bridge edges are single edges in a biconnected component
  compind <- which(sapply(bicc$component_edges, length)==1)
  if (exclude_degree_1) {
    degr <- igraph::degree(igr)
    ## if two-node biconnected components contain a node w/degree=1,
    ## exclude from list of bridge edges
    ## TODO: all vs any
    subind <- sapply(bicc$components[compind], function(v) {all(degr[v]!=1)})
    if (length(subind) == 0) return(matrix())
    else subind <- which(subind)
    compind <- compind[subind]
  }

  bridgeind <- unlist(lapply(bicc$component_edges[compind], as.numeric))
  igraph::get.edgelist(igr)[bridgeind,,drop=FALSE]
}

complementEmat <- function(A, subemat) {
  emat <- adj2emat(A)
  emat[!(emat[,1] %in% subemat[,1] & emat[,2] %in% subemat[,2]),]
}

adj2emat <- function(A, undirected=TRUE) {
  diag(A) <- 0
  emat <- arrayInd(which(A!=0), dim(A))
  if (undirected) {
    emat <- t(apply(emat, 1, sort))
  }
  unique(emat)
}

enforceEdges2 <- function(graph, e, r=1, ...) {
  ## find bridge edges to protect:
  currE <- sum(graph)/2
#  if (currE == e) return(graph)
  if (currE < e) {
    if (isTRUE(attr(graph, 'graph')=="scale_free")) {
      ## attach hubs that aren't already neighbors
      deg <- rowSums(graph)
      ## query top n hub nodes
      hubs <- head(order(deg, decreasing=TRUE), min(5+r, nrow(graph)/10))
      # which hubs are not already connected
      igr <- adj2igraph(graph)
      hubneighs <- lapply(hubs, function(h) as.numeric(igraph::neighbors(igr, h)))
      ## other hubs to consider
      newneighs <- lapply(1:length(hubneighs), function(i) hubs[!(hubs %in% c(hubneighs[[i]], hubs[i])) ])
      names(newneighs) <- hubs
      newneighs <- Filter(length, newneighs)
      if (length(newneighs)==0) {
        return(SpiecEasi:::enforceE(graph, e, ...))
      }
#      print('foo')
      node1 <- as.numeric(names(newneighs)[1])
      if (length(newneighs[[1]]) ==1) {
        node2 <- newneighs[[1]]
      } else {
        node2 <- sample(newneighs[[1]], 1, prob=deg[newneighs[[1]]])
      }
      graph[node1, node2] <- graph[node2, node1] <- 1
      return(enforceEdges2(graph, e, r=r+1, ...))
    }

    return(SpiecEasi:::enforceE(graph, e, ...))
  }

  maxIt <- 100
  it <- 1
  while (it <= maxIt && currE != e) {
    bridgE <- bridge_edges(graph, TRUE)
    if (isTRUE(nrow(bridgE) > 0) && !any(is.na(bridgE))) {
      nonbridgE <- complementEmat(graph, bridgE)
    } else {
      nonbridgE <- adj2emat(graph)
    }

#    if (isTRUE(nrow(nonbridgE) == 0)) break

    ## remove a single edge at random
    to_rm <- nonbridgE[sample(1:nrow(nonbridgE), 1),]
    graph[to_rm[1], to_rm[2]] <- graph[to_rm[2], to_rm[1]] <- 0
    # if (ncol(connected.components(graph))>2) {
    #   print(to_rm)
    #   graph[to_rm[1], to_rm[2]] <- graph[to_rm[2], to_rm[1]] <- 1
    #   return(list(graph, bridgE, nonbridgE, to_rm))
    # }
    currE <- sum(graph)/2
    it <- it+1
  }
  return(graph)
}

enforceEdges <- function(graph, e, ...) {
  if (isTRUE(attr(graph, 'graph')=="cluster")) {
    graph <- tryCatch({
      cc <- connected.components(graph)
      subgraphs <- lapply(1:ncol(cc), function(i) graph[cc[,i],cc[,i], drop=FALSE])
      subE <- sapply(subgraphs, sum)/2
      eInd <- order(subE)
      totalE <- sum(subE)
      nHubs <- length(subgraphs)
      targetE <- rmultinom(1, size=e, prob=rep(1/nHubs, nHubs))[,1]
      teInd <- order(targetE)
      for (i in 1:nHubs) {
        ind <- which(cc[,eInd[i]])
        graph[ind, ind] <- SpiecEasi::enforceE(graph[ind, ind], targetE[teInd[i]])
      }
      graph
     }(graph, e),
     error=function(err) SpiecEasi:::enforceE(graph, e, ...)
    )

 } else {
    graph <- SpiecEasi:::enforceE(graph, e, ...)
 }
  return(graph)
}

enforceDegree <- function(graph, maxdeg, r=1, maxr=maxdeg+1) {
  if (isTRUE(attr(graph, 'graph')=="cluster")) {
    cc <- connected.components(graph)
    subgraphs <- lapply(1:ncol(cc), function(i) graph[cc[,i],cc[,i], drop=FALSE])
    submaxdegs  <- sapply(subgraphs, function(x) max(rowSums(x)))
    subvsize <- sapply(subgraphs, nrow)
    sg_ind   <- order(submaxdegs, subvsize, decreasing=TRUE)[1]
    ind <- c(which(cc[,sg_ind]), which(rowSums(graph)==0))
    graph[ind, ind] <- enforceDegree(graph[ind, ind], maxdeg, r, maxr)
    deg <- rowSums(graph)
    if (max(deg) < maxdeg) {
      ## try whole graph in 'normal' mode
      attr(graph, 'graph') <- NULL
      graph <- enforceDegree(graph, maxdeg, r, maxr)
      attr(graph, 'graph') <- 'cluster'
    }
    return(graph)
  }

  deg <- rowSums(graph)
  if (max(deg) == maxdeg || r > maxr) {
    if (r > maxr) message('could not satisfy maximum degree')
    return (graph)
  } else if (max(deg) > maxdeg) {
    hub <- which.max(deg)
    neighbors <- which(graph[hub,] != 0)
    neigh     <- neighbors[sample(which(deg[neighbors] == min(deg[neighbors])), 1)]
    graph[hub, neigh] <- graph[neigh, hub] <- 0
    deg[c(hub,neigh)] <- Inf
    deg[deg==0] <- Inf
    deg[deg==maxdeg] <- Inf
    newneigh <- sample(which(!is.infinite(deg)), 1)
    graph[newneigh, neigh] <- graph[neigh, newneigh] <- 1
    return(enforceDegree(graph, maxdeg, r+1, maxr))
  } else {
    if (isTRUE(attr(graph, 'graph')=="band")) {
      ## make sure hub is close to the 'end'
      hubs <- which(deg==max(deg))
      hub <- hubs[which.min(pmin(length(deg) - hubs, hubs))]
    } else {
      ## pick any hub
      hub <- which.max(deg)
    }
    notneighbors <- setdiff(which(graph[hub,] == 0), hub)
    ## max degree is unreachable, all neighbors are exhausted
    if (length(notneighbors)==0) return(graph)
#    deg[deg==0] <- Inf

    if (isTRUE(attr(graph, 'graph')=="band")) {
      ## prefer isolated nodes as neighbors
      minnodes <- which(deg[notneighbors] == min(deg[notneighbors]))
      hubdist <- abs(notneighbors[minnodes]-hub)
      nind <- tryCatch(sample(minnodes, 1, prob=1/hubdist), error=function(e) sample(minnodes, 1))
    } else {
      ## pick any minimal node as the new neighbor
      nind <- sample(which(deg[notneighbors] == min(deg[notneighbors])), 1)
    }
    newneigh <- notneighbors[nind]

    # deg[deg==maxdeg] <- Inf
    # deg[c(hub,newneigh)] <- Inf

    ## protect bridge edges
    bridgE <- bridge_edges(graph, TRUE)
    if (isTRUE(nrow(bridgE) > 0) && !any(is.na(bridgE))) {
      nonbridgE <- complementEmat(graph, bridgE)
    } else {
      nonbridgE <- adj2emat(graph)
    }
    ## protect hub edges
    emat <- nonbridgE[!(nonbridgE[,1]==hub |
                        nonbridgE[,2]==hub),,drop=FALSE]
    ## max degree is unreachable, no edges left to steal
    if (isFALSE(nrow(emat)>0)) {
      return(graph)
    }
    oldneighs <- emat[sample.int(nrow(emat), 1),]
    graph[hub, newneigh] <- graph[newneigh, hub] <- 1
    graph[oldneighs[1], oldneighs[2]] <- graph[oldneighs[2], oldneighs[1]] <- 0
    return (enforceDegree(graph, maxdeg, r+1, maxr))
  }

}

BandE <- function(p, e, dir='over') {
  edges <- cumsum(p - 1:(p-1))
  diff  <- edges - e
  if (dir == 'over')
    edges[which(diff >= 0)[1]]
  else
    edges[tail(which(diff <= 0), 1)]
}

getNet <- function(type, p, e, d, dens=.04, ...) {
    if (type == "maxdeg") {
      q <- 1.28 #(log(100)-log(2))/(2*log(log(100)))
      maxdeg <- floor(log(p)^q)
      g <- igraph::sample_degseq(degseq <- rand_vect(p, 2*e, maxdeg), method='simple.')
      graph <- as.matrix(igraph::get.adjacency(g))
      class(graph) <- "graph"
    } else if (type == 'band') {
      d <- round(d*.9)
      e <- round(dens*p*(p-1)/2)
      supe <- BandE(p, e, 'over')
      graph <- make_graph(type, p, supe, enforce=FALSE)
      graph <- enforceEdges2(graph, e)
      graph <- enforceDegree(graph, d)
    } else if (type == 'cluster') {
      d <- round(d)
      e <- round(dens*p*(p-1)/2)
      hubs  <- max(2, round(e/p))
      graph <- make_graph(type, p, max(e, 2.5*p),
                          numHubs=hubs, enforce=FALSE)
      graph <- enforceEdges2(graph, e)
      graph <- enforceDegree(graph, d)
    } else if (type == 'scale_free') {
      d <- round(d*1.1)
      e <- round(dens*p*(p-1)/2)
      graph <- make_graph(type, p, p+1, enforce=FALSE)
      graph <- enforceEdges2(graph, e)
      graph <- enforceDegree(graph, d)
    }

    G <- diag(p) - 1/p
    Prec  <- graph2prec(graph, targetCondition=2e1)
#    Prec  <- graph2prec2(graph, targetCondition=2e1)
    Cov   <- prec2cov(Prec)

    spread <- svd(Prec)$d[1]/norm(Prec, 'M')
    d <- max(rowSums(abs(sign(Prec))))
    eta <- spread/d
    S     <- G%*%Cov%*%G
    Gammainv <- G%*%Prec%*%G
    L <- Prec - Gammainv
    L.svd <- svd(L)
    r <- sum(L.svd$d > 1e-6)
    incoh <- norm(M <- L.svd$u[,1:r] %*% t(L.svd$v[,1:r]), 'M')
    incoh2 <- get_incoh2(Prec)
    mu <- (p^2*incoh^2)/r
    return(list(graph=graph, Prec=Prec, Cov=Cov, S=S, d=d, maxdeg=NULL,
                spread=spread, d=d, eta=eta, incoh=incoh, incoh2=incoh2, mu=mu))
}

get_incoh2 <- function(Prec) {
  p <- ncol(Prec)
  U <- rowMeans(Prec) / norm(rowMeans(Prec), '2')
  UUt <- U%*%t(U)
  P <- rowMeans(diag(p) - UUt)
  P <- P/norm(P, '2')
  PPt <- P %*% t(P)
  norm(UUt+PPt, 'M')
}

make_sparse_graphical_model <- function(p=15, e=p, type="band", d=max(2, round(p*(1/4))), ...) {
  new_model(name = "sgm",
            label = sprintf("p = %s, e = %s, degree = %s, method = %s", p, e, d, type),
            params = list(p=p, e=e, d=d, type=type),
            simulate = function(p, e, type, d, nsim) {
              replicate(nsim, getNet(type, p, e, d), simplify=FALSE)
            }
  )
}

# p <- 10
# J <- matrix(1/p, p, p)
# G <- diag(p) - J
# tmp <- getNet('band', p, p-1, 2)
# #tmp$Prec <- diag(p) #rbind(cbind(c(1, 1,rep(0,p-2)), diag(0, nrow=p-1, ncol=p-1)), 0) # #1.1*diag(p) # matrix(1, p,p ) #
# #tmp$Prec[1,p] <- tmp$Prec[p,1] <- 1.1
# ##tmp$Prec[p-1,2] <- tmp$Prec[2, p-1] <- 5
# ##tmp$Prec[p-2,3] <- tmp$Prec[3, p-2] <- -1
# ##tmp$Cov <- solve(tmp$Prec)
#
# #tmp$graph <- abs(sign(tmp$Prec)) - diag(p)
# d <- max(rowSums(abs(sign(tmp$Prec))))
#
# etamin <- 1/p
# etamid <- 1/d
# etamax <- 1
#
# spread <- svd(tmp$Prec)$d[1]/norm(tmp$Prec, 'M')
# eta <- spread/d
# tmp$Cov   <- prec2cov(tmp$Prec)
# S     <- G%*%tmp$Cov%*%G
# Gammainv <- G%*%tmp$Prec%*%G
# L <- tmp$Prec - Gammainv
# L.svd <- svd(L)
# r <- sum(L.svd$d > 1e-6)
# incoh <- norm(M <- L.svd$u[,1:r] %*% t(L.svd$v[,1:r]), 'M')
# mu <- (p^2*incoh^2)/r
#
# #print("sparse")
#
# ##print(spread <= etamin*d)
# ##print(spread <= etamid*d)
# ##print(spread <= etamax*d)
# #eta <- spread/d
# #print(eta > etamax || eta < etamin)
# #print("lowrank")
#
# #mumax <- p/r
# #mumin <- 1
#
# #print(mu > mumax || mu < mumin)
# JPPJ <- J%*%tmp$Prec + tmp$Prec%*%J
# JPPJ.svd <- svd(JPPJ)
# incoh.max <- norm(M <- ((U <- JPPJ.svd$u[,1:r,drop=FALSE]) %*% t(V <- JPPJ.svd$v[,1:r,drop=FALSE])), 'M')*sqrt(r)
#
# Pj <- tmp$Prec%*%j
# Up <- Pj / norm(Pj, '2')
# ones <- rep(1, p)
# P <- (diag(p) - Up%*%t(Up))%*%ones
# #P <- P/norm(P, '2')
#
# Vp <- ones / sqrt(p)
# Q <- (diag(p) - Vp%*%t(Vp))%*%Pj
#
#
#
# incoh.max2 <- max(rowSums(cbind(Up, P)^2))
# incoh.max3 <- max(rowSums(U^2))
#
# ##JPJ <- -J%*%tmp$Prec%*%J
# ##JPJ.svd <- svd(JPJ)
# ##a <- JPJ.svd$u[,1] * sqrt(JPJ.svd$d[1])
# ##b <- JPJ.svd$v[,1] * sqrt(JPJ.svd$d[1])
# ##Sd <- diag(JPPJ.svd$d[1:r], r, r)
# ##K <- Sd + t(U)%*%a %*%t(t(V)%*%b)
# ##K.svd <- svd(K)
# ##Up <- K.svd$u
# ##Vp <- K.svd$v
# ##kappa <- 1 - (rot <- t(Up)%*%Vp)[1,1]
#
#
# #alpha <- 2*sqrt(mu*r*eta*d/p)


# graph2prec2 <- function(graph,  thetaLims=c(2,3),
#                         targetCondition=100, epsBin=1e-2,
#                         numBinSearch=100) {
#
#   Theta <- graph #+ diag(ncol(graph))
#   attributes(Theta)$graph <- NULL
#   maxDeg <- max(rowSums(Theta))
#
#   for (i in 1:maxDeg) {
#     ind <- which(i == rowSums(Theta))
#     for (j in ind) {
#       Theta[j,Theta[,j]!=0] <- sample(c(1,-1)[(1:i %% 2)+1])
#     }
#   }
#
#   Tiu <- SpiecEasi:::triu(Theta, TRUE)
#   Tiu <- runif(length(Tiu), thetaLims[1], thetaLims[2])*Tiu
#   Theta <- triu2diag(Tiu, 1)
#   eigVals <- eigen(Theta)$values
#   minEig  <- min(eigVals)
#   maxEig  <- max(eigVals)
#   n <- ncol(Theta)
#   if (minEig < 1e-2) Theta <- Theta + abs(minEig)*diag(n)
#   diagConst <- SpiecEasi:::.binSearchCond(Theta, targetCondition, numBinSearch, epsBin)
#   Theta <- Theta + diagConst*diag(n)
#   return(Theta)
# }
