load('QMPphyfilt.RData')

phy@tax_table@.Data[grepl("[a-z]__$", phy@tax_table@.Data)] <- NA
tmptax <- apply(phy@tax_table@.Data, 1, function(x) tail(x[!is.na(x)],1))
tmptax <- gsub("\\[|\\]", "", stringr::str_match(tmptax, "^[a-z]__(.*)$")[,2])
tmptax[duplicated(tmptax)] <- paste(tmptax[duplicated(tmptax)], '2')
taxa_names(phy) <- tmptax

## filter
phy.f <- prune_taxa(rowSums(sign(phy@otu_table@.Data))>0, phy)
otus <- t(phy.f@otu_table@.Data)
otus.f <- t(apply(otus,1,norm_to_total))

sample_dirmnom <- function(x, size) {
  n <- nrow(x)
  t(sapply(1:n, function(i) extraDistr::rdirmnom(1, size[i], x[i,])))
}

otus.f <- t(apply(otus+1/10, 1, SpiecEasi::norm_to_total))
qmp <- pmax(round(otus.f * phy@sam_data@.Data[[27]]), 1)

get_qmp_subset_std <- function(p, index) {
  Sig <- cov(log(qmp[,index]))
  G <- diag(p) - 1/p
  Prec <- MASS::ginv(Sig)
  spread <- svd(Prec)$d[1]/norm(Prec, 'M')
  d <- max(rowSums(abs(sign(Prec))))
  eta <- spread/d
  Cov   <- prec2cov(Prec)
  S     <- G%*%Cov%*%G
  Gammainv <- G%*%Prec%*%G
  U <- rowMeans(Prec) / norm(rowMeans(Prec), '2')
  UUt <- U%*%t(U)
  P <- rowMeans(diag(p) - UUt)
  P <- P/norm(P, '2')
  PPt <- P %*% t(P)
  incoh2 <- norm(UUt+PPt, 'M')
  L <- Prec - Gammainv
  L.svd <- svd(L)
  r <- sum(L.svd$d > 1e-6)
  incoh <- norm(M <- L.svd$u[,1:r] %*% t(L.svd$v[,1:r]), 'M')
  mu <- (p^2*incoh^2)/r

  list(index=index, otus=otus[,index], qmp=qmp[,index],
       depths=rowSums(otus), Sig=Sig,
       eta=eta, spread=spread, incoh2=incoh2, mu=mu)
}

get_qmp_subset_pnorm <- function(p, index, eps=1e-7) {
  otus.f <- t(apply(otus, 1, SpiecEasi::norm_to_total))
  otus.f[otus.f==0] <- eps
  otus.f <- t(apply(otus.f, 1, SpiecEasi::norm_to_total))
  qmp <- pmax(round(otus.f * phy@sam_data@.Data[[27]]), 1)
  list(index=index, otus=otus[,index], qmp=qmp[,index],
       depths=rowSums(otus), Sig=cov(log(qmp[,index])))
}

get_qmp_subset_dmult <- function(p, index, eps=1/10) {
  qmp <- sample_dirmnom(otus+eps, phy@sam_data@.Data[[27]])
  list(index=index, otus=otus[,index], qmp=qmp[,index],
       depths=rowSums(otus), Sig=cov(log(qmp[,index]+eps)))
}

gmean <- function(x) exp(mean(log(x[x!=0])))

make_qmp_subset <- function(p, sim_type) {
  new_model(name = "QMP-subset",
            label = sprintf("p = %s; method = %s;",
                            p, sim_type),
            params = list(p=p, type=sim_type),
            simulate = function(p, type, nsim) {
              lapply(1:nsim, function(i) {
                  ind <- sample(1:ncol(otus), p,
                          prob=norm_to_total(colMeans(otus>0)))
                  switch(type,
                    std  =get_qmp_subset_std(p, ind),
                    pnorm=get_qmp_subset_pnorm(p, ind),
                    dmult=get_qmp_subset_dmult(p, ind)
                  )
              })
            }
  )
}
