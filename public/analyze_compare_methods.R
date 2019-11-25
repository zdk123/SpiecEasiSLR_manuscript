#!/bin/env Rscript

library(dplyr)
source('analyze_net_funs.R')

currdir <- getwd()

projects <- sapply(list.dirs(recursive=FALSE),
              list.dirs, recursive=FALSE) %>%
              unlist %>% unname

d_list <- lapply(projects, function(p) {
  setwd(p)
  print(p)
  if (!file.exists('net_dist.RData')) {
    tmp <- tryCatch( {
      source('../../compare_nets.R')
    }, error=function(e) NULL, warning=function(w) NULL)
    if (is.null(tmp)) out <- NULL
  } else {
    tmp <- 1
    load('net_dist.RData')
  }
  setwd(currdir)
  if (!is.null(tmp)) {
    list(edge_dist=as.matrix(edge_dli),
         edge_dense=as.matrix(edge_dense),
         edge_ppos=as.matrix(edge_ppos),
         degree_dist=as.matrix(degree_dli),
         btw_dist=as.matrix(btw_dli)) -> out
  } else out <- NULL
  out
})

names(d_list) <- projects
d_list <- d_list[!sapply(d_list, is.null)]

d_df <- reshape2::melt(d_list)

lsort <- function(x, y) x[order(ordered(x, levels = llevs))]

llevs <- c("se.poi", "se.is", "se.gl", "se.mb","se.slr")
d_df[,1:2] <- data.frame(t(apply(d_df[,1:2], 1, lsort, y= c(llevs,'1'))), stringsAsFactors=FALSE)
d_df <- unique(d_df)
d_df <- d_df[d_df[,1] != d_df[,2],]
d_df$Var2[d_df$Var2=='1'] <- d_df$Var1[d_df$Var2=='1']

d_df$Var1 <- factor(d_df$Var1, levels=llevs)
d_df$Var2 <- factor(d_df$Var2, levels=llevs)
d_df$Var3 <- 'x'
library(ggplot2)
library(ggridges)

tmp <- d_df %>% dplyr::filter(L2 %in% c('edge_dist'))
p1 <- ggplot(aes(x=L2,y=value), data=tmp) +
        facet_grid(Var2~Var1,drop=FALSE) +
        stat_summary(fun.y=median, geom='point', shape='_', size=8, color='black', alpha=.5) +
        geom_jitter(size=.3, width=.1, color='red', alpha=.6) +

        scale_fill_brewer(palette='Set1') +
        theme_bw() +
        theme(axis.text.y = element_text(size=6),
              axis.text.x = element_text(size=6)) +
        guides(fill=FALSE)
      #  coord_flip()


tmp <- d_df %>% dplyr::filter(L2 %in% c('degree_dist', 'btw_dist'))
p2 <- ggplot(aes(y=L2,x=(1-value)/2, fill=L2), data=tmp) +
      #  theme_ridges() +
        facet_grid(Var1~Var2,drop=FALSE) +
        stat_density_ridges(alpha=.6, size=.2, bandwidth=.06,
            panel_scaling=FALSE, jittered_points=TRUE, point_shape=16,
            point_alpha=.6, point_size=.3,
            position = position_points_sina(rel_min=.1, rel_max=.65),
            quantile_lines=TRUE, quantiles=2) +
        scale_fill_brewer(palette='Set2') +
        theme_bw() +
        theme(axis.text.y = element_text(size=6 , angle = 90, hjust = 1),
              axis.text.x = element_text(size=6, angle=45, hjust=1)) +
        guides(fill=FALSE) #+
      #  coord_flip()



tmp <- d_df %>% dplyr::filter(L2 %in% c('edge_dense'))
p3 <- ggplot(aes(x=value, fill=L2), data=tmp) +
        facet_grid(Var3~Var2,drop=FALSE) +
        scale_fill_manual(values='black') +
        geom_histogram(bins=12, alpha=.5) +
        theme_bw() +
        theme(axis.text.x = element_text(size=6, angle=45, hjust=1),
              axis.text.y = element_text(size=5)) +
        guides(fill=FALSE)


tmp <- d_df %>% dplyr::filter(L2 %in% c('edge_ppos'))
p4 <- ggplot(aes(x=value, fill=L2), data=tmp) +
        facet_grid(Var3~Var2,drop=FALSE) +
        scale_fill_manual(values='black') +
        geom_histogram(bins=12, alpha=.5) +
        theme_bw() +
        theme(axis.text.x = element_text(size=6, angle=45, hjust=1),
              axis.text.y = element_text(size=5)) +
        guides(fill=FALSE)



png('dev/null')
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)

#dev.off()
maxW <- grid::unit.pmax(g1$widths[2:5], g2$widths[2:5], g3$widths[2:5], g4$widths[2:5])

g1$widths[2:5] <- as.list(maxW)
g2$widths[2:5] <- as.list(maxW)
g3$widths[2:5] <- as.list(maxW)
g4$widths[2:5] <- as.list(maxW)

pdf('../figures/net_dists.pdf', height=8, width=3.5)
gridExtra::grid.arrange(g1,g2,g3,g4, heights=c(3,3,1.2,1.2))
dev.off()
