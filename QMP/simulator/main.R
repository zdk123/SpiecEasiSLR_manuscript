library(simulator)
library(SpiecEasi)
library(pulsar)
library(Matrix)
library(phyloseq)
library(ggplot2)
source("model_functions.R")
source('method_functions.R')
source('eval_functions.R')


ps <- round(rev(pulsar::getLamPath(205, 16, 7, log=TRUE)))
sim_type <- list('dmult')

psock <- list(socket_names=4,
              libraries=loadedNamespaces())

sim <- new_simulation('qmp-spiec-easi', "QMP SPIEC-EASI") %>%
  generate_model(make_qmp_subset, p = as.list(ps),
                 sim_type=sim_type,
                 vary_along = c("p", "sim_type"),
                 seed=10010) %>%
  simulate_from_model(nsim = 1, index=1:10, parallel=psock) %>%
  run_method(list(se.slr,se.gl,coat,crude), parallel=psock) %>%
  evaluate(list(minFNorm,  oracleNorm,  SigNorm,
                minFNorm2, oracleNorm2, SigNorm2,
                minFNorm3, oracleNorm3, SigNorm3))

sim <- load_simulation('qmp-spiec-easi')

evdf   <- as.data.frame(evals(sim))
evdf$p <- as.numeric(stringr::str_match(as.character(evdf$Model), "p_(.*)/")[,2])
evdf$eps <- as.numeric(stringr::str_match(as.character(evdf$Model), "eps_([\\.0-9]*)/")[,2])
evdf$mod2 <- basename(as.character(evdf$Model))
cols <- c(segl='#66c2a5', semb='#fc8d62', crude='#e78ac3',
          sepois='#a6d854', seslr='#8da0cb', coat='#e5c494',
          seis='#ffd92f')


plot_fun <- function(Var) {
  ggplot(aes(x=p, color=Method), data=evdf) +
        # facet_wrap(~p, scales='free', nrow=1) +
       # geom_boxplot(aes_string(y=Var), outlier.color='grey60') +
        stat_summary(aes_string(y=Var), fun.y = mean,
                     geom = "line", size=1.5, alpha=1) +
        stat_summary(aes_string(y=Var), fun.data = mean_cl_normal,
                     geom = "errorbar", width=.075) +
        scale_x_log10(breaks=ps) +
        scale_color_manual(values=cols) +
        theme_bw(base_size=16)
}

p1 <- plot_fun('Fnorm/SigNorm')
p2 <- plot_fun('Fnorm_shrink/SigNorm_shrink')
p3 <- plot_fun('OracleFnorm')
p4 <- plot_fun('OracleFnorm_shrink')

pdf('../../figures/sim_qmp.pdf', width=4.5, height=5)
egg::ggarrange(p1 + ylab("Relative Fro Norm"),
               p3,ncol=1,newpage=FALSE)
dev.off()
