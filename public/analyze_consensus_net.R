## check pairwise consistency between nets ##
source("scripts/analyze_net_funs.R")
library(dplyr)
library(Matrix)
library(ggplot2)

projects <- grep('/', list.files(pattern='[^/]NetFits.RData', recursive=TRUE), value=TRUE)

# nets_li <- lapply(projects, load_nets, sign=TRUE)
# save(nets_li, file="signed_nets_list.RData")
load('signed_nets_list.RData')
names(nets_li) <- dirname(projects)
nets_li <-  purrr::transpose(nets_li)


gr2edges <- function(x) {
  edges <- Matrix::summary(x)
  edges[,1] <- rownames(x)[edges[,1]]
  edges[,2] <- rownames(x)[edges[,2]]
  edges[,1:2] <- t(apply(edges[,1:2], 1, sort))
  edges <- as.data.frame(edges)
  edges$edge <- paste(edges[,1], edges[,2], sep=".")
  tmp <- split(edges, edges$x)
  names(tmp) <- ifelse(names(tmp)=="1", "pos", "neg")
  tmp
}


edges_li <- lapply(nets_li, function(nets) purrr::transpose(lapply(nets, gr2edges)))

edges2adj <- function(edges) {
  pos_edges <- table(unlist(lapply(edges$pos, function(x) x$edge)))
  neg_edges <- table(unlist(lapply(edges$neg, function(x) x$edge)))

  ## build consensus weighted adjancency matrix
  taxa <- sort(unique(unname(
        c(unlist(lapply(edges$pos, function(x) x[,1:2])),
          unlist(lapply(edges$neg, function(x) x[,1:2])))
        )))

  rebuild_adj <- function(x) {
    if (is.null(x)) {
      return(Matrix::Matrix(data=0, ncol=length(taxa), nrow=length(taxa)))
    }
    sign(sparseMatrix(i=match(x$i, taxa),
                 j=match(x$j, taxa),
                 x=x$x, dims=c(length(taxa), length(taxa))))
  }

  edges_full <- list(pos=lapply(edges$pos, rebuild_adj),
                     neg=lapply(edges$neg, rebuild_adj))

  pos_ <- Reduce(`+`, x=edges_full$pos)/
            length(edges_full$pos)
  neg_ <- Reduce(`+`, x=edges_full$neg)/
            length(edges_full$neg)

  adj <- forceSymmetric(pos_ + neg_)
  rownames(adj) <- colnames(adj) <- taxa
  adj
}
## TODO: get consensus nets for all methods
adj_li  <- lapply(edges_li, edges2adj)
adj_fli <- lapply(adj_li, function(adj) {
                 adj[abs(adj) < 5/26] <- 0
                 deg <- Matrix::rowSums(adj)
                 adj[deg!=0, deg!=0]
            })


subtax <- unique(unlist(lapply(adj_fli, rownames)))

cols <- c('ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
library(data.table)
gg_tax_all <- fread('gg_13_5_taxonomy_split.txt.gz', col.names=cols)
gg_tax_all[, ID := as.character(ID)]
setkey(gg_tax_all, ID)
gg_tax_all[gg_tax_all==""] <- NA
gg_tax <- gg_tax_all[subtax]
setkey(gg_tax, ID)


## color key @ order level:
vcol <- c(
  Clostridiales='#df89ff',
  Bacteroidales='#73c000',
  Enterobacteriales='#00c4ff',
  Lactobacillales='#4c463e',
  Bifidobacteriales='#ff8805',
  Actinomycetales='#ff5584',
  Erysipelotrichales='#00bd94',
  Burkholderiales='#d3b3b0',
  Pasteurellales='#fbff38',
  default='#c0c0c0')

library(igraph)
consadj2igraph <- function(adj) {
  tmpadj <- forceSymmetric(adj)
  gg_sub <- gg_tax[rownames(tmpadj), -(1:2)]
  SpiecEasi::adj2igraph(abs(tmpadj), #(tmpadj+1)/2 * abs(tmpadj),
            rmEmptyNodes=FALSE,
            vertex.attr=as.pairlist(gg_sub),
            edge.attr=list(sign=Matrix::summary(sign(t(tmpadj)))[,3])
            )
}
igr_li <- lapply(adj_fli, consadj2igraph)

largestCC <- function(igr) {
  ## get the largest connected component
  cc <- components(igr)
  vind <- which(cc$membership== which.max(cc$csize))
  induced_subgraph(igr, vind)
}
igrCC_li <- lapply(igr_li, largestCC)

rgr <- rgexf::igraph.to.gexf(igr_li[['se.slr']])
print(rgr, file='net_consens.gexf')


getMetrics <- function(igr_li) {

if (inherits(igr_li, 'igraph')) {
    igr_li <- list(igr_li)
}

# assortmat <- sapply(igr_li, function(igr) {
#   vert_phy <- as.data.frame(gg_tax[V(igr)$name,'Phylum'])
#   vert_fam <- as.data.frame(gg_tax[V(igr)$name,'Family'])
#   vert_phy[vert_phy==""] <- NA
#   vert_fam[vert_fam==""] <- NA
#   c(
#   phylum=1-assortativity_nominal(induced_subgraph(igr, which(vert_phy!="")),
#                        as.factor(vert_phy[!is.na(vert_phy)])),
#   order=1-assortativity_nominal(induced_subgraph(igr, which(vert_fam!="")),
#                        as.factor(vert_fam[!is.na(vert_fam)]))
#   # phylum=assortativity_nominal(igr,
#   #                      as.factor(vert_phy)),
#   # order=assortativity_nominal(igr,
#   #                      as.factor(vert_fam))
#   )
# })
#
# centrmat <- sapply(igr_li, function(igr) {
#   c(
#   degree=centr_degree(igr)$centralization,
#   betw=centr_betw(igr)$centralization
#   )
# })
#
# countmat <- sapply(igr_li, function(igr) {
#   c(
#     vcount=vcount(igr),
#     ecount=ecount(igr)
#   )
# })
#
# rbind(
#   cbind(reshape2::melt(centrmat), metric='centralization'),
#   cbind(reshape2::melt(assortmat), metric='assortativity'),
#   cbind(reshape2::melt(countmat), metric='count')
# )
# }
#
# netmet_df <- getMetrics(igr_li)

## get all module subgraphs for se.slr
igr_li_mod <- list(se.slr=lapply(unique(modindex), function(i) induced_subgraph(igr_li$se.slr, which(modindex==i))))

## get order
assort <- function(igr) assortativity_nominal(igr, as.factor(V(igr)$Phylum))



netmetCC_df <- getMetrics(igrCC_li)


cols <- c(otu='#BEBEBE', se.mb='#fdbf6f', se.gl='#66c2a5', se.slr='#8da0cb', se.is='#cab2d6', se.poi='#b15928')

p1 <- ggplot(aes(Var1, value, fill=Var2, group=Var2), data=netmet_df) +
  facet_wrap(~metric, scales='free') +
  geom_bar(stat='identity', position='dodge', color='black') +
  scale_fill_manual(values=cols)


p2 <- ggplot(aes(Var1, value, fill=Var2, group=Var2), data=netmetCC_df) +
  facet_wrap(~metric, scales='free') +
  geom_bar(stat='identity', position='dodge', color='black') +
  scale_fill_manual(values=cols)


# egg::ggarrange(p1,p2, nrow=2)

## decompose graph into triangles
method <- 'se.slr'
tri <- cliques(igr_li[[method]], min=3, max=3)
itri <- lapply(tri, induced_subgraph, graph=igr_li[[method]])
itax <- lapply(itri, function(x) {
    y <- gg_tax[V(x)$name,-(1:2)] %>%
         arrange(Phylum, Class, Order, Family, Genus)
    rownames(y) <- V(x)$name
    y
})

## find triangles with a negative edge
neg_ind <- sapply(itri, function(gr) any(E(gr)$sign==-1))


as.data.frame(t(sapply(itax, function(x) x$Family))) %>%
    mutate(has_neg=neg_ind) %>%
    group_by(V3,V2,V1) %>%
    summarize(n=n(), has_neg=sum(has_neg)) %>%
    na.exclude %>% ungroup %>%
    mutate(n_uniq=apply(.[,1:3], 1, function(x) length(unique(x)))) %>%
    arrange(-n_uniq, -n) -> tritax


## find triangles with different families
ind <- sapply(itri, function(x) length(unique(V(x)$Family)) > 2)
# itax[ind]


set.seed(10010)
modindex <- membership(cluster_fast_greedy(igr_li$se.slr, weights=abs(E(igr_li$se.slr)$sign)))


## Write edge list for SLR graph
slr_edgeli <- summary(forceSymmetric(adj_fli$se.slr))
slr_edgeli[,1] <- rownames(adj_fli$se.slr)[slr_edgeli[,1]]
slr_edgeli[,2] <- rownames(adj_fli$se.slr)[slr_edgeli[,2]]
slr_edgeli[,4] <- sign(slr_edgeli[,3])
slr_edgeli[,3] <- abs(slr_edgeli[,3]*length(projects))
colnames(slr_edgeli) <- c("gg1", "gg2", "ndatasets", "sign")
slr_edgeli <- slr_edgeli[order(slr_edgeli$ndatasets, decreasing=TRUE),]

slr_edgeli$module1 <- modindex[slr_edgeli$gg1]
slr_edgeli$module2 <- modindex[slr_edgeli$gg2]
slr_edgeli$taxonomy1 <- apply(gg_tax[slr_edgeli$gg1,-1],1,paste, collapse=";")
slr_edgeli$taxonomy2 <- apply(gg_tax[slr_edgeli$gg2,-1],1,paste, collapse=";")
## add modules
write.table(slr_edgeli, "slr_consensedglist.txt", sep="\t", row.names=FALSE, quote=FALSE)



set.seed(10010)
plot(igr_li$se.slr, vertex.size=5,
     vertex.color=vcol[V(igr_li$se.slr)$Order],
     vertex.label.cex=.5,
     edge.color=as.factor(E(igr_li$se.slr)$sign),
     vertex.label=modindex)
## subnet module index 6
#induced_subgraph(igr_li$se.slr, V(igr_li$se.slr)$name[modindex==6])
#induced_subgraph(igr_li$se.slr, do.call('union',ego(igr_li$se.mb, nodes=n))$name)

get_datasets <- function(vnames, method='se.slr') {
  emat <- combn(vnames,2)
  n <- length(vnames)

  lapply(edges_li[[method]], function(signednet) {
    edges <- apply(cbind(emat, emat[2:1,]), 2, paste, collapse=".")
    tmp <- sapply(1:n-1, function(i)
      which(sapply(signednet, function(x) {
        any(edges[c(1,n+1)+i] %in% x[,4])
    })))
    names(tmp) <- edges[1:n]
    tmp
  })
}
#neg_itax <- lapply(tri[neg_ind], function(x) gg_tax[x$name,])


##igrsub <- induced_subgraph(igr_li[[method]], which(V(igr_li[[method]])$name %in% unique(row.names(do.call('rbind', neg_itax)))))



### plot all graphs
plot_net <- function(adj, rank='Order', metanodes, ...) {
  ind <- which(rowSums(abs(adj))!=0)
  adj <- as.matrix(adj[ind,ind])
  tax <- gg_tax_all[rownames(adj)]$Order
  col <- vcol[tax]
  col[is.na(col)] <- vcol['default']
  eCol <- ifelse(sign(adj)>0,'darkgreen','firebrick')
  eLty <- matrix(1, nrow(adj), ncol(adj))
  eCex <- matrix(1, nrow(adj), ncol(adj))
  eLwd <- matrix(1, nrow(adj), ncol(adj))
  if (!missing(metanodes)) {
    eCol[metanodes,] <- 'grey70'
    eCol[,metanodes] <- 'grey70'
    eLty[metanodes,] <- 3
    eLty[,metanodes] <- 3
    eLwd[metanodes,] <- 1/3
    eLwd[,metanodes] <- 1/3
  }
  plot(network::network(adj), usearrows=FALSE,
       edge.col=eCol, edge.lty=eLty,
       vertex.col=col, ...)
}

igr_alli <- lapply(nets_li$se.slr, function(adj) {
  tmpadj <- forceSymmetric(adj)
  gg_sub <- as.data.frame(gg_tax_all[rownames(tmpadj),-(1:2)])

  mod6 <- rownames(adj) %in% names(which(modindex==6))

  SpiecEasi::adj2igraph(abs(tmpadj),
    vertex.attr=c(list(name=rownames(tmpadj),
              mod6=rownames(adj) %in% names(which(modindex==6))),
                    as.pairlist(gg_sub)),
    edge.attr=list(sign=Matrix::summary(sign(t(tmpadj)))[,3])
    )
})

igrsub_alli <- lapply(igr_alli, largestCC)

## write GGMP
rgr <- rgexf::igraph.to.gexf(igrsub_alli[[3]])
print(rgr, file=file.path(dirname(projects[3]), 'Subnet.gexf'))


### layout
# library(ForceAtlas2)
# layout_fa2 <- function(igr) {
#   tmpig <- igr
#   E(tmpig)$weight <- ifelse(E(tmpig)$sign==1,1,.5)
#   layout.forceatlas2(tmpig, directed=FALSE, iterations=1000, gravity=1, delta=.5, plotstep=500)
# }
#
# layout_list <- lapply(igrsub_alli, layout_fa2)

pdf('../figures/slr_nets.pdf')
for (i in 1:length(nets_li$se.slr)) {
  labs <- V(igrsub_alli[[i]])$name
  tmp <- as.matrix(get.adjacency(igrsub_alli[[i]], attr='sign'))
  colnames(tmp) <- rownames(tmp) <- labs
  plot_net(tmp, edge.lwd=.4)
}
dev.off()


whiten <- function(x, d=3) {
  xsvd <- svd(x)
  xsvd$v[,1:d] %*% diag(xsvd$d[1:d]) %*% t(xsvd$v[,1:d])
}

## mod6 resid matrix
icov2pcor <- function(S) {
  -cov2cor(S) + 2*diag(ncol(S))
}

source('seriate_slrresids.R')
pdf("../figures/slr_resids_mod6.pdf", width=3, height=3)
for (i in 1:length(igr_alli)) {
  ind <- os[[i]]$heir[[1]]$order
#  ind <- as.vector(os[[i]]$pca[[1]])
  subind <- which(V(igr_alli[[i]])$mod6[ind])

  if (length(subind) > 0) {
    labs <- V(igr_alli[[i]])$name[ind][subind]
    L <- resids[[i]]$Lcor[ind,ind][subind,subind]
    plot(gimage(icov2pcor(L), labs, main=basename(names(resids)[i])))
  }
}
dev.off()

## plot mod6 subindex
pdf('../figures/slr_nets_mod6.pdf', height=5, width=10)
for (i in 1:length(igr_alli)) {
  ind <- os[[i]]$heir[[1]]$order
  subind <- which(V(igr_alli[[i]])$mod6[ind])

  if (length(subind) > 0) {
    labs <- V(igr_alli[[i]])$name[ind][subind]
    S <- as.matrix(resids[[i]]$S[ind,ind][subind,subind])
    colnames(S) <- rownames(S) <- labs
    L <- as.matrix(resids[[i]]$L[ind,ind][subind,subind])
    R <- S-L
    colnames(R) <- rownames(R) <- labs
    R[abs(R)<=.001] <- 0
    par(mfrow=c(1,2))
    plot_net(sign(icov2pcor(S)), vertex.cex=4,
             main=basename(names(igr_alli)[i]),
             label=labs, label.cex=.6, edge.lwd=3,
             label.pad=0) -> layout
    plot_net(sign(icov2pcor(R)), vertex.cex=4, coord=layout,
      main=basename(names(igr_alli)[i]), edge.lwd=3,
      label=labs, label.cex=.6, label.pad=0)
  }
}
dev.off()


## get GGMP correlates
i <- 'qiita/GGMP'
ind <- os[[i]]$heir[[1]]$order
subind <- which(V(igr_alli[[i]])$mod6[ind])
L <- as.matrix(resids[[i]]$L[ind,ind][subind,subind])
S <- as.matrix(resids[[i]]$S[ind,ind][subind,subind])
load('qiita/GGMP/rob_pca.RData')
load('qiita/GGMP/GGMPphyfilt.RData')
X <- t(SpiecEasi::clr(t(phy_taxfilt@otu_table@.Data+1), 1))

Lsvd <- svd(L)
dind <- Lsvd$d>3e-3
loadings <- diag(sqrt(1/Lsvd$d[dind])) %*% t(Lsvd$v[,dind])
scores <- X[,ind][,subind] %*% t(loadings)
colnames(scores) <- 1:4

source('analyze_rob_pca_funs.R')

col_fact <- read.table('qiita/GGMP/map_fact_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]
col_num <- read.table('qiita/GGMP/map_num_cols.txt', header=FALSE, stringsAsFactors=FALSE)[,1]

xcorr <- cross_corr(scores, map, col_num, col_fact) %>%
    arrange(desc(R^2)) %>% group_by(Var1) %>% top_n(3, wt=R^2)


scor <- cor(scores, X[,ind][,subind]) #%>% melt(varnames=c("PC", "OTU"), value.name='R')
scor[abs(scor)<=.5] <- 0

Saug <- rbind(cbind(sign(solve(S))*abs(sign(S)), t(scor)), cbind(diag(4), scor))
Saug <- Saug[22:1,22:1]
diag(Saug) <- 0

coord <- plot_net(Saug, metanodes=1:4, vertex.cex=3)
coord[1:4,2] <- predict(lm(cy~cx, data=as.data.frame(coord[1:4,])))

pdf('qiita/GGMP/mod6_net.pdf')
plot_net(Saug, metanodes=1:4, vertex.cex=4, coord=coord)
dev.off()
