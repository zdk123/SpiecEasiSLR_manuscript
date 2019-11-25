  ## check pairwise consistency between nets ##
source("scripts/analyze_net_funs.R")
projects <- grep('/', list.files(pattern='[^/]NetFits.RData', recursive=TRUE), value=TRUE)

library(dplyr)

if (file.exists("all_nets_list.RData")) {
  load('all_nets_list.RData')
} else {
  nets_li <- lapply(projects, load_nets)
  names(nets_li) <- dirname(projects)
  save(nets_li, file='all_nets_list.RData')
}

nets_li <-  purrr::transpose(nets_li)

ham_dist <- lapply(nets_li, function(nets) {
  proxy::dist(nets, method=X_dist, f=edge_dist)
})

deg_dist <- lapply(nets_li, function(nets) {
  proxy::dist(nets, method=X_dist, f=degree_cordist)
})

otu_dist <- proxy::dist(nets_li$se.slr, method=node_dist)

dist2mat <- function(x) {
  as(forceSymmetric(as.matrix(x)), 'dsCMatrix')
}

## distance plots
dist2mat2 <- function(x) {
  summary(dist2mat(x))
}

dist2df <- function(ham_dist) {
  ham_dfli <- lapply(ham_dist, dist2mat2)
  ham_df <- reshape2::melt(ham_dfli, id.var=c('i','j'), meas.var='x')
  colnames(ham_df) <- c('i', 'j', 'X', 'distance', 'method')
  ham_df$project <-  paste(projects[ham_df$i],  projects[ham_df$j], sep="-")
  ham_df$method <- factor(ham_df$method, level=c('se.poi','se.is','se.mb','se.gl','se.slr'))
  ham_df
}

library(ggplot2)
library(ggridges)
library(dplyr)

cols <- c(otu='#BEBEBE', se.mb='#fdbf6f', se.gl='#66c2a5', se.slr='#8da0cb', se.is='#cab2d6', se.poi='#b15928')

ham_df <- dist2df(ham_dist)
ham_df %>% group_by(project) %>%
           dplyr::filter(distance==min(distance)) -> ham_min
tmp <- as.data.frame(t(do.call('cbind', strsplit(ham_min$project, "-"))))
colnames(tmp) <- c('p1', 'p2')
rownames(tmp) <- NULL
ham_min <- bind_cols(ham_min, tmp)

p1 <- ggplot(aes(x=distance), data=ham_df) +
          geom_point(aes(y=method), data=ham_min, size=2, shape='|') +
          geom_density_ridges(aes(y=method, fill=method, height=..density..), alpha=.8, trim = TRUE, stat = "density", bw=0.034) +
          stat_density_ridges(aes(y=method), quantile_lines = TRUE, quantiles = 2, alpha=0, bandwidth=0.034) +
          scale_fill_manual(values=cols) +
          scale_y_discrete(expand=c(0,0)) +
          xlim(0,1) +
          guides(fill="none") +
          theme_ridges()


bind_rows(ham_min, tibble(method='dummy', distance=.5)) -> ham_min2
ham_min2$method <- factor(ham_min2$method, levels=c(levels(ham_min$method), 'dummy'))
p2 <- ggplot(aes(x=as.numeric(method), y=distance, size=method!='dummy'), data=ham_min2) +
        stat_summary(geom='linerange', fun.data=function(x) data.frame(y=0, ymin=0, ymax=length(x)), size=5) +
        scale_x_discrete(expand=c(0,0)) +
        scale_size_manual(c('TRUE'=.1, 'FALSE'=5)) +
        coord_flip() + theme_void()

p1 <- egg::ggarrange(p1,p2,nrow=1, widths=c(4,1), draw=FALSE)
pdf('../figures/dataset_dist.pdf', width=5, height=4)
print(p1)
dev.off()


(ham_df %>% group_by(method) %>% summarize(med=median(distance)) -> tmp)
