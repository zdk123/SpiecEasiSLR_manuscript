library(simulator)
library(SpiecEasi)
library(pulsar)
library(Matrix)
library(phyloseq)
library(ggplot2)
r <- 43
source('process.R')
## fit true nets or load existing RData object
#source('fit_truenets.R')
load('TrueNets.RData')


source("model_functions.R")
source('method_functions.R')
source('eval_functions.R')


nmax <- nsamples(phy_train)
ns <- round(exp(seq(log(nmax/60), log(nmax), length.out=6)))

psock <- list(socket_names=10,
	      libraries=loadedNamespaces())

# sim <- new_simulation('agp-spiec-easi', "AGP SPIEC-EASI") %>%
#  generate_model(make_AGP_sub, n = as.list(ns),
#                 vary_along = "n", seed=10010) %>%
#  simulate_from_model(nsim = 1, index=1:5) %>%
#  run_method(list(se.is, se.poi, se.gl, se.mb, se.slr)) %>%
#   evaluate(list(maxf1, maxf1_2, hamming, hamming_2, zeronorm, zeronorm_2))
sim <- load_simulation('agp-spiec-easi') #%>%

evdf   <- as.data.frame(evals(sim))
evdf$n <- stringr::str_match(evdf$Model, "n_([0-9]+)")[,2] %>%
              as.numeric

cols <- c('se-gl'='#66c2a5', 'se-mb'='#fc8d62', crude='#e78ac3', 'se-poi'='#a6d854', 'se-slr'='#8da0cb', coat='#e5c494', 'se-is'='#ffd92f')

plot_fun <- function(var) {
  ggplot(aes_string(x='n', y=var, color='Method'), data=evdf) +
        # geom_boxplot(aes_string(y=Var), outlier.color='grey60') +
        stat_summary(fun.y = mean, geom = "line", size=2, alpha=.7) +
        stat_summary(fun.data = mean_sdl, geom = "errorbar", width=.1) +
        scale_x_log10(breaks=ns) +
        # scale_y_continuous(limits=c(.45,1)) +
        scale_color_manual(values=cols) +
        theme_bw(base_size=12)
}


p1 <- plot_fun('hamming/(Znorm*2)')
p2 <- plot_fun('F1')

pdf('../../../../figures/AGP_cons.pdf')
egg::ggarrange(p1 + ylab("Relative Hamming"),
               p2 + ylab("Max F1"),
               ncol=1, newpage=FALSE)
dev.off()
