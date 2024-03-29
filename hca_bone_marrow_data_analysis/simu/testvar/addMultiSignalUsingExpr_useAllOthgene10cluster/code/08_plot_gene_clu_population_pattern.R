# /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/code/08_plot_gene_clu_population_pattern.R
rm(list=ls())
library(ggplot2)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
library(parallel)
library(splines)
library(viridis)
library(here)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

# rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/'
# pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/clu/'

rdir <- here('hca', 'simu', 'testvar', 'addMultiSignalUsingExpr', 'result')
pdir <- here('hca', 'simu', 'testvar', 'addMultiSignalUsingExpr', 'plot', 'clu')

Res <- readRDS(paste0(rdir, '/EM/1.rds'))
pt = Res$pseudotime
clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/clu/geneclu.rds')


# ---------------------------------
# fit sample level (not preferable)
# ---------------------------------
pred <- Res$predict.values[names(sort(Res$fdr[Res$fdr < 0.05])), names(pt)]
agg <- t(sapply(unique(clu), function(i) colMeans(pred[clu == i, ])))
rownames(agg) <- paste0('cluster', unique(clu))

library(ggplot2)
pd <- reshape2::melt(agg)
pd$pseudotime = pt[pd[,2]]
colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
pd$cluster <- factor(as.character(pd$cluster), levels =  paste0('cluster', 1:max(clu)))

pdf(paste0(pdir, 'gene_cluster_pattern.pdf'), width = 4.5, height = 3.2)  
ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
  geom_smooth() +
  theme_classic() +
  scale_color_manual(values =  brewer.pal(10,'Set3')) +
  xlab('Pseudotime') +
  ylab('Averaged Fitted Centered Expression')
dev.off()

# ------------------------------------------
# plot each sample seperately: inappropriate
# ------------------------------------------
cellanno = Res$cellanno
rownames(cellanno) = cellanno[,1]
cellanno = cellanno[colnames(agg), ]

design = Res$design

agg.s = agg[, cellanno[,2] %in% rownames(design[design[,2] ==0, ])]
pd <- reshape2::melt(agg.s)
pd$pseudotime = pt[as.character(pd[,2])]
colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
# pd$cluster <- factor(as.character(pd$cluster), levels =  paste0('cluster', 1:max(clu)))
pdf(paste0(pdir, 'gene_cluster_expr_along_pseudotime_all_samples.pdf'), width = 4.5, height = 3.2)  
ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
  geom_smooth() +
  theme_classic() +
  scale_color_manual(values =  brewer.pal(10,'Set3')) +
  xlab('Pseudotime') +
  ylab('Averaged Fitted Centered Expression')
dev.off()
# --------------------------------------------------
# fit population pattern for the cell-level clusters
# --------------------------------------------------
fit <- lapply(1:max(clu), function(i){
  ag = names(clu[clu==i])
  fit <- lapply(ag, function(g){
    get_population_fit(Res, variable = 'group', gene = g)
  })
  
  fit1 <- colMeans(t(sapply(fit, function(g){
    g[[1]]
  })) )
  fit0 <- colMeans(t(sapply(fit, function(g){
    g[[2]]
  })))
  rbind(data.frame(group = 1, pseudotime = seq(1, length(fit1)), expr=fit1, cluster = i, stringsAsFactors = FALSE),data.frame(group = 0, pseudotime = seq(1, length(fit1)), expr=fit0, cluster = i, stringsAsFactors = FALSE))
})

pd <- do.call(rbind, fit)
saveRDS(pd, paste0(pdir, 'clu_population_pattern_pd.rds'))

pd$group = as.factor(pd$group)
library(ggplot2)
pdf(paste0(pdir, 'clu_population_pattern_pd.pdf'), width = 8, height = 8)
ggplot() + geom_point(data = pd, aes(x = pseudotime, y = expr, color = group, group = group)) +
  theme_classic() +
  facet_wrap(~cluster, scale = 'free') + 
  ylab('Cluster averaged population fit on centered values')
dev.off()


# ---------------------------------------------------
# fit population pattern for clustering and then plot
# ----------------------------------------------------
## fit population level curves for differential genes on variable "group"
fit <- t(sapply(names(Res$fdr[Res$fdr < 0.05]), function(g){
  tmp <- get_population_fit(Res, variable = 'group', gene = g)
  c(tmp[[1]], tmp[[2]])
}))
## clustering
set.seed(12345)
clu <- kmeans(fit,10,iter.max = 1000)$cluster
clu <- sort(clu)
## average the population fit for genes within a same cluster
fit1 <- sapply(1:max(clu), function(i){
  colMeans(fit[names(clu[clu==i]), 1:(ncol(fit)/2),drop = FALSE])
})
fit0 <- sapply(1:max(clu), function(i){
  colMeans(fit[names(clu[clu==i]), (ncol(fit)/2 + 1):ncol(fit),drop = FALSE])
})
rownames(fit1) <- rownames(fit0) <- paste0(1:nrow(fit1))
colnames(fit1) <- colnames(fit0) <- paste0('cluster', 1:ncol(fit1))
## make plot data
pd1 <- cbind(group = 1, reshape2::melt(fit1))
pd0 <- cbind(group = 0, reshape2::melt(fit0))
pd <- rbind(pd1, pd0)
colnames(pd) <- c('group', 'pseudotime', 'cluster', 'expr')
saveRDS(pd, paste0(pdir, '/clu_based_on_population_pattern_pd.rds'))

pd$group = as.factor(pd$group)
library(ggplot2)
pdf(paste0(pdir, '/clu_based_on_population_pattern_pd.pdf'), width = 8, height = 6)
ggplot() + geom_point(data = pd, aes(x = pseudotime, y = expr, color = group, group = group)) +
  theme_classic() +
  facet_wrap(~cluster, scale = 'free') + 
  ylab('Cluster averaged population fit on centered values')
dev.off()


