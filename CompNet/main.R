library(simulator) # this file was created under simulator version 0.2.0
#library(SpiecEasi)
devtools::load_all('~/NYU/SpiecEasi')
library(latex2exp)
library(ggplot2)
library(dplyr)
library(Matrix)
source("simulator/model_functions.R")
source("simulator/method_functions.R")
source("simulator/eval_functions.R")

## init

name_of_simulation <- "lr-spiec-easi"
set.seed(10010)
## @knitr main
#densities <- as.list(rev(pulsar::getLamPath(.1, .0034, 3, log=TRUE)))
ps <- round(rev(pulsar::getLamPath(150, 16, 7, log=TRUE)))
pargs <- list(socket_names = 5, libraries = c("huge", 'SpiecEasi', 'Matrix'))

types <- c("scale_free", "cluster", "band")
sim <- new_simulation('lr-spiec-easi', "Low Rank SPIEC-EASI") %>%
  generate_model(make_sparse_graphical_model, p = as.list(ps),
                 type=as.list(types), vary_along = c("p", "type")) %>%
 simulate_from_model(nsim = 1, index=1:10, parallel = pargs) %>%
  run_method(list(coat, se.slr, se.gl), parallel = pargs) %>%
# sim <- load_simulation('lr-spiec-easi') %>%
   evaluate(list(hamming, maxf1, spread, incoh, incoh2, perfectf1, spread.incoh, zeronorm.icov, zeronorm.cov, varrowmeans.icov, maxdeg, mu, eta, oracle_size))

##sim <- load_simulation('lr-spiec-easi')


### plots
evdf <- as.data.frame(evals(sim))
evdf$p <- as.numeric(gsub("p_", "", sapply(strsplit(as.character(evdf$Model), "\\/"), '[[', 2)))
evdf$model <- gsub("type_", "", sapply(strsplit(as.character(evdf$Model), "\\/"), '[[', 3))

df2 <- reshape2::melt(evdf, measure.vars=c("Zero_cov", "Zero_icov", "maxdeg"))
evdf$incoh_ub <- sqrt(1/evdf$p)
df3 <- reshape2::melt(evdf, measure.vars=c("incoh", "incoh2", "spread"))

p2 <- ggplot(aes(x=p, y=value, color=variable), data=df3) +
       facet_wrap(~model) +
       stat_summary(fun.y = mean, geom = "line", size=1.3) +
       stat_summary(fun.data = mean_sdl, geom = "errorbar", width=.075) +
       scale_x_log10(breaks=ps[seq(2,9,2)]) +
       scale_y_log10() +
       ylab("") +
       theme_bw(base_size=14)

p3 <- ggplot(aes(x=p, y=value, color=variable), data=df2) +
       facet_wrap(~model) +
       stat_summary(fun.y = mean, geom = "line", size=1.3) +
       stat_summary(fun.data = mean_sdl, geom = "errorbar", width=.075) +
       scale_x_log10(breaks=ps[seq(2,9,2)]) +
       scale_y_log10(breaks=10^(1:4),
            labels = scales::trans_format("log10", scales::math_format(10^.x ) )) +
       theme_bw(base_size=14) +
       ylab("") +
       scale_color_brewer(palette='Set1',
          labels=list(
              TeX("$\\frac{1}{2} || \\Omega || _{0}$"),
              TeX("$\\frac{1}{2} || \\Omega^{-1}|| _{0}$"),
              TeX("$d_{max}$")))

cols <- c(segl='#66c2a5', seslr='#8da0cb', coat='#e5c494')
p4 <- ggplot(aes(x=p, y=F1, color=Method), data=evdf) +
       facet_wrap(~model) +
       stat_summary(fun.y = mean, geom = "line", size=1.3) +
       stat_summary(fun.data = mean_se, geom = "errorbar", width=.075) +
       scale_x_log10(breaks=ps) +
       scale_color_manual(values=cols) +
       theme_bw(base_size=14)

p5 <- ggplot(aes(x=p, y=hamming/Zero_icov, color=Method), data=evdf) +
      facet_wrap(~model) +
      stat_summary(fun.y = mean, geom = "line", size=1.3) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width=.075) +
      scale_x_log10(breaks=ps) +
      scale_color_manual(values=cols) +
      theme_bw(base_size=14)


pdf('../figures/sim_compnet.pdf', height=9, width=8)
egg::ggarrange(p4+ylab("max F1"),
               p5+ylab("relative Hamming"), p2,p3, ncol=1, newpage=FALSE)
par(mfrow=c(1,3))
for (i in 0:2) {
 plot(network::network(draws(sim)[[5+(i*8)]]@draws[[1]]$graph), usearrows=FALSE, displayisolates=FALSE, vertex.cex=2)
}
dev.off()
