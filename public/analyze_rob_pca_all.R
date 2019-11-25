library(dplyr)

currdir <- getwd()

projects <- sapply(list.dirs(recursive=FALSE),
              list.dirs, recursive=FALSE) %>%
              unlist %>% unname

robpca_list <- lapply(projects, function(p) {
  setwd(p)
  print(p)
  if (!file.exists('rob_pca.RData')) {
    tmp <- tryCatch( {
      source('scripts/analyze_rob_pca_funs.R')
      suppressWarnings(source('robPCA.R'))
      save(pvals, keep, scores, map, xCors, dims=dims, file='rob_pca.RData')
      1
    }, error=function(e) NULL)
    if (is.null(tmp)) out <- NULL
  } else {
    tmp <- 1
    load('rob_pca.RData')
  }
  setwd(currdir)
  if (!is.null(tmp)) {
    list(pvals=pvals,
         keep=keep,
         xcors=xCors,
         dims=dims,
         scores=scores,
         map=map) -> out
  } else out <- NULL
  out
})
names(robpca_list) <- projects

## Collect some statistics
num_features <- sapply(robpca_list, function(x) {
  x$pvals$n_features
})

num_keep <- sapply(robpca_list, function(x) {
  length(x$keep[[1]])
})

num_samples <- sapply(robpca_list, function(x) {
  x$dims[1]
})

num_otus <- sapply(robpca_list, function(x) {
  x$dims[2]
})

unexplained_comp <- sapply(robpca_list, function(x) {
  x$keep$unexplained_components
})
n_comp <- sapply(robpca_list, function(x) {
  x$keep$n_comp
})


sq <- function(x)  if(is.null(x)) NULL else x^2
## effect size distributions ##
eff <- lapply(robpca_list,
          function(x) do.call('rbind', lapply(x$pvals[[1]],
            function(y) {
              tmp <- rbind(y$quali[,1,drop=FALSE],
                        sq(y$quanti[,1,drop=FALSE]))
             data.frame(R2=unname(tmp), variable=rownames(tmp))
        })))
## get the max effect size for each dataset/variable for all components
effmax <- lapply(eff, function(x) x %>% group_by(variable) %>% summarize(R2=max(R2)))

catmap <- lapply(basename(projects), function(n) {
  tmp <- xlsx::read.xlsx('significant_features.xlsx', sheetName=n, header=FALSE)
  colnames(tmp) <- c('variable', 'category')
  tmp[,1:2]
})
names(catmap) <- basename(projects)

eff <- lapply(1:length(eff), function(i) {
  suppressWarnings(left_join(eff[[i]], catmap[[i]], by='variable'))
})
names(eff) <- names(catmap)

effmax <- lapply(1:length(effmax), function(i) {
  suppressWarnings(left_join(effmax[[i]], catmap[[i]], by='variable'))
})
names(effmax) <- names(catmap)

library(ggplot2)
eff_tbl <- reshape2::melt(effmax)
colnames(eff_tbl)[3:4] <- c('R2', 'dataset')

amax <- function(x) data.frame(maxc=x[which.max(abs(x))], i=which.max(abs(x)))
## hist of max correlations
corr_df <- do.call('rbind', lapply(1:length(robpca_list), function(i) data.frame(robpca_list[[i]]$xcors, dataset=names(robpca_list)[i]))) %>%
  rename(variable=L1, i=Var1) %>%
  mutate(qval=p.adjust(pval, method="bonferroni"),
           R2=R^2)

left_join(corr_df,
    data.frame(dataset=names(n_comp),
               n_comp=unname(n_comp))) %>%
left_join(
    data.frame(dataset=names(num_features),
               n_feat=unname(num_features))) -> corr_df

cat_map <- read.delim('data/feature_cats.txt')
dplyr::left_join(corr_df, cat_map) %>%
  mutate(category=ifelse(is.na(category), 'other', as.character(category))) -> corr_df

catnames <- sort(names(table(corr_df$category)))
catnames <- c(catnames[5], catnames[1:4])
corr_df$category <- factor(corr_df$category, levels=catnames)
colmap <- setNames(c('grey70', RColorBrewer::brewer.pal(5, 'Set1')), c(catnames))


corr_df %>% group_by(dataset, variable) %>% arrange(desc(R2)) %>% top_n(1, wt=R2) -> corrsub_df

corr_df %>% group_by(dataset, i) %>% arrange(desc(R2)) %>% top_n(1, wt=R2) -> corrsub2_df

corrsub_df %>% dplyr::filter(qval<=10^-5) %>%
    group_by(dataset, n_comp, n_feat) %>%
    summarize(n=n_distinct(i)) -> tmp

p1 <- lattice::xyplot(tmp$n_comp~tmp$n/tmp$n_feat, type=c('p','r'),
                      xlab="Significant Features",
                      ylab="Latent Components")


p5 <- ggplot(aes(R2, group=category, fill=category), data=eff_tbl) +
        geom_histogram(alpha=.7, bins=25) +
        scale_y_continuous(min=eff_tbl$R2) +
        theme_bw() + scale_fill_manual(values=colmap)

p6 <- ggplot(aes(R2, group=category, fill=category), data=corr_df) +
        facet_wrap(~dataset, scales='free_y') +
        geom_histogram(alpha=.7, bins=8) +
        scale_fill_manual(values=colmap) +
        theme_bw()

library(ggridges)
p7 <- ggplot(aes(R2, y=0, fill=category, height=..count..), data=corrsub_df) +
         facet_grid(category~., scales='free_y') +
         geom_density_ridges(stat='binline', binwidth=0.0875) +
        scale_fill_manual(values=colmap, guide=FALSE) +
        scale_x_continuous(breaks=seq(0,1,.2)) +
        coord_cartesian(xlim=c(0.04,.96)) +
        theme_ridges(font_size=12) +
        theme(panel.spacing = grid::unit(-.2, "lines"))

pdf('tmp1.pdf', height=4, width=4)
print(p1)
#gridExtra::grid.arrange(p1,p2,ncol=2)
print(p5)
print(p7)
dev.off()

pdf('tmp2.pdf', height=7, width=12)
gridExtra::grid.arrange(p6, ncol=1)
dev.off()

system('qpdf --empty --pages tmp* -- ../figures/robpca.pdf')
system('rm tmp*')
## variable categories
eff_tbl %>% group_by(category, variable) %>%
      summarize(n=n(), R2=mean(R2)) %>%
      arrange(desc(n), desc(R2)) %>%
      dplyr::filter(n>=3) -> eff_top
print(xtable::xtable(eff_top), file='../figures/robpca_tab.tex')



### plot some of the top correlations
corrsub_df %>% dplyr::filter(category!='other') %>%
  group_by(category) %>% arrange(category, -R2) %>%
  top_n(4, wt=R2) -> topcorr_df

get_data <- function(x) {
  x <- as.data.frame(t(x))
  dataset <- as.character(x$dataset)
  variable <- as.character(x$variable)
  scores <- robpca_list[[dataset]]$scores
  map <- robpca_list[[dataset]]$map
  data.frame(score=scores[,x$i],
           feature=map[,variable])
}

topscores <- apply(topcorr_df, 1, get_data)

## plot all
plt_li <- vector('list', nrow(topcorr_df))
pdf('../figures/robpca_top.pdf', height=3, width=5)
for (i in 1:nrow(topcorr_df)) {
  print(as.data.frame(topcorr_df[i,]))
  tmpdf <- topscores[[i]]

  if (class(tmpdf$feature) == 'factor') {
    tmpdf$feature <- factor(tmpdf$feature,
                            levels=tmpdf %>% group_by(feature) %>%
                               summarize(med=mean(score)) %>%
                               arrange(-med) %>% .$feature)


    plot(plt_li[[i]] <- ggplot(aes(feature, score), data=tmpdf) +
      geom_point(position=position_jitter(width=.15, seed=10010),
                 size=.6, alpha=.5) +
      stat_summary(fun.data=mean_cl_normal,
                   geom='errorbar', width=.3, alpha=.5) +
      stat_summary(fun.y=mean, geom='point', color='red',
                   size=3, alpha=.6, shape="|") +
      xlab("") +
      ylab(sprintf("PC%s", topcorr_df[i,'i'])) +
      theme_bw() + coord_flip() +
      ggtitle(paste(topcorr_df[i,]$category,
                    basename(topcorr_df[i,]$dataset),
                    topcorr_df[i,]$variable, sep=" - "))
    )
  } else {
    plot(plt_li[[i]] <- ggplot(aes(feature, score), data=tmpdf) +
         geom_point(size=.5, alpha=.6) +
          xlab("") +
          ylab(sprintf("PC%s", topcorr_df[i,'i'])) +
          theme_bw() + coord_flip() +
          ggtitle(paste(topcorr_df[i,]$category,
                        basename(topcorr_df[i,]$dataset),
                        topcorr_df[i,]$variable, sep=" - "))
        )
  }
}
dev.off()

## aligned subset
ind <- c(2, 7, 9, 13)
plts <- egg::ggarrange(plots=plt_li[ind], ncol=1, draw=FALSE, newpage=FALSE)
ggsave('../figures/robpca_top2.pdf', plts, height=9, width=4)
