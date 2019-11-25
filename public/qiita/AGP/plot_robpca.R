## Example Viz
library(ggplot2)
library(dplyr)
load('scores.RData')
load('rob_pca.RData')

cats <- c('zero_tax', 'depth', 'center_project_name', 'vegetable_frequency', 'library_construction_protocol', 'instrument_model', 'country_residence', 'census_region', 'economic_region')


#ind <- apply(abs(xCors), 2, function(x) order(x, decreasing=TRUE)[1:2])
xCors %>% mutate(R2=R) %>%
          filter(L1%in%cats) %>%
          arrange(desc(R2)) %>%
          group_by(L1) %>%
          top_n(2, wt=R2) %>%
          select(L1, Var1) %>%
          tidyr::spread(L1, Var1) %>% as.data.frame -> ind

mean_sd <- function(x) {
  mean <- mean(x)
  sd <- sqrt(var(x))
  data.frame(y = mean, ymin = mean - sd, ymax = mean + sd)
}

## Zero tax / model
i <- 1
scores_df <- data.frame(X1=scores[,ind[,cats[i]]], zero_tax=map[,cats[i]], model=map[,cats[6]])

# cor test
mod <- cor.test(scores_df$X1, scores_df$zero_tax)
mod_df <- data.frame(x=Inf,y=-Inf,
          hjustvar=1.1, vjustvar=-.1,
          label=paste(
            sprintf("r^2 = %g", signif(mod$estimate^2,2)),
            sprintf("pval %s", format.pval(mod$p.value,4)),
          sep="\n"))

p1 <- ggplot(aes(X1, zero_tax), data=scores_df) +
#  stat_bin_hex(size=.06, alpha=.5) +
  geom_point(aes(color=model), size=.06, alpha=.2) +
  xlab(sprintf("PC%s", ind[,cats[i]][1])) +
  ylab("Zero counts") +
  geom_text(aes(x,y,label=label,hjust=hjustvar,vjust=vjustvar), data=mod_df, size=3) +
  theme_bw() +
  ggtitle(cats[i])


## Depth
i <- 2
scores_df <- data.frame(X1=scores[,ind[,cats[i]]], zero_tax=map[,cats[i]], model=map[,cats[6]])

# cor test
mod <- cor.test(scores_df$X1, scores_df$zero_tax)
mod_df <- data.frame(x=Inf,y=-Inf,
          hjustvar=1.1, vjustvar=-.1,
          label=paste(
            sprintf("r^2 = %g", signif(mod$estimate^2,2)),
            sprintf("pval %s", format.pval(mod$p.value,4)),
          sep="\n"))

p15 <- ggplot(aes(X1, zero_tax), data=scores_df) +
#  stat_bin_hex(size=.06, alpha=.5) +
  geom_point(aes(color=model), size=.06, alpha=.2) +
  xlab(sprintf("PC%s", ind[,cats[i]][1])) +
  ylab("Depth") +
  geom_text(aes(x,y,label=label,hjust=hjustvar,vjust=vjustvar), data=mod_df, size=3) +
  theme_bw() +
  ggtitle(cats[i])


## Veg freq
i <- 4
scores_df <- data.frame(X1=scores[,ind[,cats[i]]], meta=map[,cats[i]], meta2=map[,cats[6]]) %>% filter(!is.na(meta))
scores_df$meta <- factor(scores_df$meta,
                         levels=scores_df %>% group_by(meta) %>%
                            summarize(med=mean(X1)) %>%
                            arrange(-med) %>% .$meta)
## recode Never
levels(scores_df$meta)[levels(scores_df$meta)=='Never'] <- "Rarely (less than once/week)"
## recode cat names for display
levels(scores_df$meta) <- gsub("\\(", "\\\n\\(", levels(scores_df$meta))

## model relationship
mod <- anova(mod0<-lm(X1~meta, data=scores_df))
mod_df <- data.frame(x=-Inf,y=Inf,
          hjustvar=1.1, vjustvar=-.1,
          label=paste(
            sprintf("r^2 = %g", signif(cor(predict(mod0), scores_df$X1)^2,2)),
            sprintf("pval %s", format.pval(mod$Pr[1],4)),
          sep="\n"))


p2 <- ggplot(aes(meta, X1), data=scores_df) +
#  facet_wrap(~meta2) +
  geom_point(position=position_jitter(width=.15, seed=10010), size=.06, alpha=.1) +
  stat_summary(fun.data=function(x) mean_cl_normal(x, mult=10), geom='errorbar', width=.5) +
  stat_summary(fun.y=mean, geom='point', color='red', size=5, alpha=.6, shape="|") +
  xlab("") +
  ylab(sprintf("PC%s", ind[,cats[i]][1])) +
  geom_text(aes(x,y,label=label,hjust=hjustvar,vjust=vjustvar), data=mod_df, size=3) +
  theme_bw() + coord_flip() +
  ggtitle(cats[i])


## Prep protocol
i <- 5
map[,cats[i]] <- as.factor(gsub(" amplification of 16SrRNA V4|\"", "", map[,cats[i]]))
scores_df <- data.frame(X1=scores[,ind[,cats[i]]], meta=map[,cats[i]], meta2=map[,cats[2]])
## recode for display
levels(scores_df$meta) <- gsub("515", "\\\n515", levels(scores_df$meta))
scores_df$meta <- factor(scores_df$meta,
                            levels=scores_df %>%
                            group_by(meta) %>%
                            summarize(med=mean(X1)) %>%
                            arrange(-med) %>% .$meta)
## model relationship
mod <- anova(mod0<-lm(X1~meta, data=scores_df))
mod_df <- data.frame(x=-Inf,y=-Inf,
          hjustvar=-.1, vjustvar=-.1,
          label=paste(
            sprintf("r^2 = %g", signif(cor(predict(mod0), scores_df$X1)^2,2)),
            sprintf("pval %s", format.pval(mod$Pr[1],4)),
          sep="\n"))


p3 <- ggplot(aes(meta, X1), data=scores_df) +
  geom_point(position=position_jitter(width=.15, seed=10010), size=.06, alpha=.1) +
  stat_summary(fun.data=function(x) mean_cl_normal(x, mult=10), geom='errorbar', width=.5) +
  stat_summary(fun.y=mean, geom='point', color='red', size=5, alpha=.6, shape="|") +
    xlab("") +
    ylab(sprintf("PC%s", ind[,cats[i]][1])) +
    geom_text(aes(x,y,label=label,hjust=hjustvar,vjust=vjustvar), data=mod_df, size=3) +
    theme_bw() + coord_flip(ylim=c(-12,12)) +
    ggtitle(cats[i])


## Center project name
i <- 3
scores_df <- data.frame(X1=scores[,ind[,cats[i]]], meta=map[,cats[i]], meta2=map[,cats[2]])
scores_df$meta <- factor(scores_df$meta,
                         levels=scores_df %>% group_by(meta) %>%
                            summarize(med=mean(X1)) %>%
                            arrange(-med) %>% .$meta)

## model relationship
mod <- anova(mod0<-lm(X1~meta, data=scores_df))
mod_df <- data.frame(x=-Inf,y=-Inf,
          hjustvar=-.1, vjustvar=-.1,
          label=paste(
            sprintf("r^2 = %g", signif(cor(predict(mod0), scores_df$X1)^2,2)),
            sprintf("pval %s", format.pval(mod$Pr[1],4)),
          sep="\n"))

p4 <- ggplot(aes(meta, X1), data=scores_df) +
  geom_point(position=position_jitter(width=.15, seed=10010), size=.06, alpha=.1) +
  stat_summary(fun.y=mean, geom='point', color='red', size=.5, alpha=.6, shape=15) +
  xlab("") +
  ylab(sprintf("PC%s", ind[,cats[i]][1])) +
  geom_text(aes(x,y,label=label,hjust=hjustvar,vjust=vjustvar), data=mod_df, size=3) +
  theme_bw() +
  theme(axis.text.y=element_text(size=5)) +
  coord_flip(ylim=c(-12,6)) +
  ggtitle(cats[i])



## Country / Region
map[,cats[7]] <- ifelse(is.na(map[,cats[8]]), as.character(map[,cats[7]]), sprintf('US-%s', as.character(map[,cats[8]])))

i <- 7
scores_df <- (data.frame(X1=scores[,ind[,cats[i]]], meta=map[,cats[i]]) %>%
  filter(!is.na(meta))->tmp) %>%
  group_by(meta) %>%
  summarize(n=n()) %>% filter(n>=200) %>%
  left_join(tmp)
scores_df$meta <- factor(scores_df$meta,
                      levels=scores_df %>% group_by(meta) %>%
                      summarize(med=mean(X1)) %>%
                      arrange(-med) %>% .$meta)
## model relationship
mod <- anova(mod0<-lm(X1~meta, data=scores_df))
mod_df <- data.frame(x=Inf,y=Inf,
          hjustvar=1.1, vjustvar=1.1,
          label=paste(
            sprintf("r^2 = %g", signif(cor(predict(mod0), scores_df$X1)^2,2)),
            sprintf("pval %s", format.pval(mod$Pr[1],4)),
          sep="\n"))

p5 <- ggplot(aes(meta, X1), data=scores_df) +
  geom_point(position=position_jitter(width=.15, seed=10010), size=.06, alpha=.1) +
  stat_summary(fun.data=function(x) mean_cl_normal(x, mult=10), geom='errorbar', width=.5) +
  stat_summary(fun.y=mean, geom='point', color='red', size=5, alpha=.6, shape="|") +
  xlab("") +
  ylab(sprintf("PC%s", ind[,cats[i]][1])) +
  geom_text(aes(x,y,label=label,hjust=hjustvar,vjust=vjustvar), data=mod_df, size=3) +
  theme_bw() +
  coord_flip(ylim=c(-13,12)) +
  ggtitle(cats[i])


pdf('../../../figures/AGP_pca.pdf', height=14, width=5)
egg::ggarrange(p1,p15,p2,p3,p4,p5, ncol=1, newpage=FALSE)
dev.off()
