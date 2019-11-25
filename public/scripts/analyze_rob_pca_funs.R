## Make robPCA from Low Rank component ##
build_pca <- function(net_rda) {
  suppressWarnings(load(net_rda))
  library(Matrix)
  rm(se.gl, se.mb, se.poi, se.is)

  ## build PCA ##
  ind <- se.slr$ebic>1e-1
  se.slr <- se.slr[ind][[which.min(se.slr$ebic[ind])]]
  L <- se.slr$est$resid[[se.slr$select$stars$opt.index]]

  X <- se.slr$est$data
  Lsvd <- svd(L)
  ind <- Lsvd$d>1e-9
  loadings <- diag(sqrt(1/Lsvd$d[ind])) %*% t(Lsvd$v[,ind])
  scores <- X %*% t(loadings)

  ## full model ##
  S <- se.slr$est$icov[[se.slr$select$stars$opt.index]]
  Rsvd <- svd(S-L)
  ind <- Rsvd$d>1e-9
  Rloadings <- diag(sqrt(1/Rsvd$d[ind])) %*% t(Rsvd$v[,ind])
  Rscores <- X %*% t(Rloadings)

  list(X=X, scores=scores, Rscores=Rscores, dims=dim(X))
}


process_map <- function(col_fact, col_num, col_date,
                rmcats=c('X.SampleID', 'Description'), rm_fact='',
                mapping_file='combined_mapping_file.txt') {
  ## get combined map ##
  map <- read.table(mapping_file, sep="\t", header=TRUE, row.names=1, quote="")

  collect <- colnames(map)
  if (!missing(col_fact))
    collect <- col_fact
  if (!missing(col_num))
    collect <- c(collect, col_num)
  if (!missing(col_date))
    collect <- c(collect, names(col_date))

  keep <- colnames(map) %in% collect
  map <- map[,keep]
  ## parse meaningful categories ##
  ncats <- apply(map, 2,  function(x) length(table(x)))
  keep <- which(ncats > 1 & !(names(ncats) %in% rmcats))

  ## Get non-redundant columns ##
  l <- lapply(map[,keep], function(X) as.numeric(factor(X, levels=unique(X))))
  m <- as.matrix(data.frame(l))
  M <- cor(m,m, use="pairwise.complete.obs")==1
  M[is.na(M)] <- FALSE
  M[lower.tri(M, diag=TRUE)] <- FALSE
  keep2 <- colnames(M)[colSums(M)==0]

  if (!missing(col_fact)) {
    col_fact <- col_fact[col_fact %in% colnames(map)]
    for (c in col_fact)
      map[,c] <- as.factor(map[,c])
  }
  if (!missing(col_num)) {
    col_num <- col_num[col_num %in% colnames(map)]
    for (c in col_num)
      map[,c] <- as.numeric(map[,c])
  }
  if (!missing(col_date)) {
    fmt <- unname(col_date)
    col_date <- names(col_date)
    fmt <- fmt[col_date %in% colnames(map)]
    col_date <- col_date[col_date %in% colnames(map)]
    for (i in 1:length(col_date)) {
        map[,col_date[i]] <- as.numeric(as.POSIXct(map[,col_date[i]], format=fmt[i]))
    }
  }

  rm_fact <- setNames(rep(NA, length(rm_fact)), rm_fact)
  rm_fact <- c(which(apply(map[,col_fact], 2, function(x) length(table(x)) >= length(x[!is.na(x)])-3)), rm_fact)
  keep3 <- keep2[!(keep2 %in% names(rm_fact))]
  map[,keep3]
}


attach_depth <- function(phy, map) {
  if (inherits(phy, 'character')) {
    load(phy)
  } else {
    phy_taxfilt <- phy
  }
  depths <- colSums(phy_taxfilt@otu_table@.Data)
  zero_tax  <- colSums(sign(phy_taxfilt@otu_table@.Data))
  map$depth <- depths
  map$zero_tax <- zero_tax
  return(map)
}


cross_corr <- function(scores, map, col_num, col_fact) {
  qual_cor <- function(x, y) {
    xna <- which(!is.na(x))
    mod <- lm(y~x)
    aov <- anova(mod)
    c(R=cor(predict(mod), y[xna]), pval=aov$Pr[1])
  }
  quant_cor <- function(x, y) {
    cormod <- cor.test(x, y,
                       use='pairwise.complete')
    c(R=unname(cormod$estimate), pval=cormod$p.value)
  }

  if (length(col_fact)>0) {
    col_fact <- col_fact[col_fact %in% colnames(map)]
    cor_fact <- parallel::mclapply(map[,col_fact], function(x) apply(scores, 2, function(y) qual_cor(x, y)), mc.cores=4)
#    cor_fact <- do.call('cbind', cor_fact)
  } else
    cor_fact <- c()

  if (length(col_num)>0) {
    col_num <- unique(c(col_num[col_num %in% colnames(map)], 'depth', 'zero_tax'))
    cor_num <- parallel::mclapply(map[,col_num], function(x) apply(scores, 2, function(y) quant_cor(x, y)), mc.cores=4)
#    cor_num <- do.call('cbind', cor_num)
  } else
    cor_num <- c()

  tmp <- c(cor_fact, cor_num)
  tmp <- reshape2::melt(lapply(Filter(length, tmp), t))
  tmp %>% dplyr::filter(Var2=='R') %>% select(-Var2) %>%
      rename(R=value) %>%
  left_join(tmp %>% dplyr::filter(Var2=='pval') %>%
      select(-Var2) %>% rename(pval=value), by=c("Var1", "L1"))
}

signif_factors <- function(scores, map, max.pval=1e-12) {
  ## Find Associations ##
  pvals <- parallel::mclapply(1:ncol(scores), function(i) condes(data.frame(scores[,i], map), num.var=1, proba=max.pval),
    mc.cores=6, mc.silent=FALSE)
  n <- ncol(scores) * ncol(map)
  list(rapply(pvals, function(x) {
    pval.adj <- p.adjust(x[,2], 'bonferroni', n)
    cbind(x[pval.adj<=max.pval,,drop=FALSE],
      adj.pval=pval.adj[pval.adj<=max.pval]
        )
    }, classes='matrix', how='replace')
  )

}


signif_features <- function(pvals) {
  quali_pval <- lapply(pvals, '[[', 'quali')
  # cat_pval   <- lapply(pvals, '[[', 'category')
  quanti_pval <- lapply(pvals, '[[', 'quanti')
  list(
    unique(unlist(sapply(c(quali_pval, quanti_pval), rownames))),
    unexplained_components = length(pvals) - (sum(  sapply(quali_pval, length) | sapply(quanti_pval, length))),
    n_comp = length(pvals)
  )
}

compare_models <- function(scores, Rscores, map) {
## compare ordination models ##
  Rord <- vegan::rda(Rscores~., data=map, na.action=na.omit)
  ord  <- vegan::rda(scores~., data=map, na.action=na.omit)

  list(
  ord=ord, Rord=Rord,
  relative_likelihood = (extractAIC(ord)[2] - extractAIC(Rord)[2])/2
  )
}

plot_pca <- function(scores, map, pvals) {
  library(tibble)
  library(dplyr)

  melt_list <- function(pli) {
    tmp <- lapply(1:length(pli), function(i) {
          tmpx <- pli[[i]]
          rownames(tmpx) <- NULL
          data.frame(
            feature=rownames(pli[[i]]),
            tmpx,
            component=rep(i, length(pli[[i]])/3 ))
          })
    Reduce(rbind, Filter(nrow, tmp))
  }

  quali_pval <- lapply(pvals, '[[', 'quali') %>%
                    melt_list %>% as_tibble
  cat_pval   <- lapply(pvals, '[[', 'category')  %>%
                    melt_list %>% as_tibble
  quanti_pval <- lapply(pvals, '[[', 'quanti') %>%
                    melt_list %>% as_tibble
###
  quali_pval %>% group_by(feature) %>% filter(order(-R2)<=3)

  ell <- function(x) {
      tryCatch(suppressMessages(as.data.frame(ellipse::ellipse(cov(x[,1:2]), centre=colMeans(x[,1:2]), lev=.3))),
      error=function(e) data.frame())
  }

  plot_qual <- function(scores, feature) {
    keep <- map[,feature] %in% names(which(table(map[,feature]) > 1))
    df <- data.frame(scores[keep,1:2], var=map[keep,feature])
    df %>% group_by(var) %>% do(ell(.)) -> ellDF
    ggplot(aes(x=X1, y=X2, z=var, group=var, color=var), data=df) +
    geom_path(aes(x=X1,y=X2, group=var), alpha=.5, data=ellDF) +
    geom_point() +
    guides(color=FALSE)
  }

}


condes <- function (donnee, num.var, weights = NULL, proba = 0.05)
{
  ## copy of FactorMineR::condes with a few memory optimizations

    cor.calc <- function(y, x, w = NULL) {
        if (is.null(w))
            w = rep(1, length(x))
        Z <- cbind(x, y)
        missing <- apply(is.na(Z), 1, any)
        Z <- Z[!missing, ]
        w <- w[!missing]
        n = sum(w)
        if (n < 3)
            n <- sum(w) * length(x)
        r = cov.wt(Z, wt = w, method = "ML", cor = TRUE)$cor[1,
            2]
        return(list(r = r, proba = pt(sqrt(n - 2) * sqrt(r^2/(1 -
            r^2)), n - 2, lower.tail = FALSE) * 2))
    }
    test.aov.w <- function(y, x, w = NULL) {
        if (is.null(w))
            w = rep(1, length(x))
        res.aov <- aov(y ~ x, weights = w, na.action = na.exclude)
        res <- summary(res.aov)[[1]]
        ddlR <- sum(w[!apply(is.na(cbind.data.frame(x, y)), 1,
            any)]) - nlevels(x)
        tabF <- c(res[1, 2]/(res[1, 2] + res[2, 2]), pf((res[1,
            3])/(res[2, 2]/(ddlR)), res[1, 1], ddlR, lower.tail = FALSE))
        Estimate <- summary.lm(res.aov)$coef[-1, 1, drop = FALSE]
        Estimate <- c(Estimate, -sum(Estimate))
        tabX <- FactoMineR:::tab.disjonctif(x)
        aux <- apply(tabX, 2, cor.calc, y, w = w)
        aux <- matrix(as.numeric(sapply(aux, unlist)), byrow = T,
            ncol = 2)
        p.value <- aux[, 2]
        resT <- cbind(Estimate, p.value)
        return(list(tabF = tabF, resT = resT))
    }
    donnee <- as.data.frame(donnee)
    is.quali <- which(!unlist(lapply(donnee, is.numeric)))
    donnee[, is.quali] <- lapply(donnee[, is.quali, drop = FALSE],
        as.factor)
    donnee <- droplevels(donnee)
    lab.sauv <- lab <- colnames(donnee)
#    quali = NULL
    if (is.null(weights))
        weights <- rep(1, nrow(donnee))
    if (sum(weights) < 3)
        weights <- weights * nrow(donnee)
    quali <- vector('numeric', length(lab))
    for (i in 1:length(lab)) {
        if (is.factor(donnee[, i])) {
            if (any(is.na(donnee[, i]))) {
                levels(donnee[, i]) <- c(levels(donnee[, i]),
                  "NA")
                donnee[, i][is.na(donnee[, i])] <- "NA"
            }
            if (levels(donnee[, i])[1] == "")
                levels(donnee[, i])[1] = "NA"
            if (i != num.var)
#                quali = c(quali, i)
                quali[i] = i
        }
    }
    quali <- quali[quali!=0]
    quanti = (1:ncol(donnee))[-c(quali, num.var)]
    if (length(quanti) == 0)
        quanti = NULL
    colnames(donnee) = lab
    result = list()
    if (!is.null(quanti)) {
        if (length(quanti) > 1) {
            tab.quanti <- apply(donnee[, quanti], 2, cor.calc,
                donnee[, num.var], w = weights)
            aux <- matrix(as.numeric(sapply(tab.quanti, unlist)),
                byrow = TRUE, ncol = 2)
        }
        else aux <- matrix(unlist(cor.calc(donnee[, quanti],
            donnee[, num.var], w = weights)), ncol = 2)
        rownames(aux) = colnames(donnee)[quanti]
        resQ = NULL
        if (NROW(aux) > 1)
            aux <- aux[rev(order(aux[, 1])), ]
        resQ <- aux[aux[, 2] < proba, , drop = FALSE]
        colnames(resQ) = c("correlation", "p.value")
        if (nrow(resQ) == 0)
            resQ = NULL
        result$quanti <- resQ
    }
    if (!is.null(quali)) {
        old.contr = options()$contrasts
        options(contrasts = c("contr.sum", "contr.sum"))
        tabF = matrix(NA, length(quali), 2)
        tabT = matrix(NA, 1, 2)
        indice.tabT = 0
        for (v in 1:length(quali)) {
            resaov <- test.aov.w(donnee[, num.var], donnee[,
                quali[v]], w = weights)
            tabF[v, ] <- resaov$tabF
            resT <- resaov$resT
            rownames(resT) <- paste(colnames(donnee)[quali[v]],
                  levels(donnee[, quali[v]]), sep=".")
#            rownames(resT) = levels(donnee[, quali[v]])
            tabT = rbind(tabT, resT)
        }
        rownames(tabF) = colnames(donnee)[quali]
        colnames(tabF) = c("R2", "p.value")
        tabT = tabT[-1, ]
        resF = resT = NULL
        if (sum(tabF[, 2] < proba) > 0)
            resF <- tabF[tabF[, 2] < proba, , drop = FALSE]
        if (!is.null(resF))
            resF <- resF[order(resF[, 2]), , drop = FALSE]
        tabT <- tabT[rev(order(sign(tabT[, 1])/tabT[, 2])), ]
        if (sum(tabT[, 2] < proba) >= 1)
            resT <- tabT[tabT[, 2] < proba, , drop = FALSE]
        result$quali = resF
        result$category = resT
        options(contrasts = old.contr)
    }
    if (is.null(result$quanti) & is.null(result$quali) & is.null(result$category))
        result = NULL
    return(result)
}
