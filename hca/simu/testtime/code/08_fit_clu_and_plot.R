# ----------------------------------------------------
# fit population pattern, clustering, and then plot
# ----------------------------------------------------
## fit population level curves for differential genes 

library(ggplot2)
library(RColorBrewer)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))
ddir <- here('hca','data','simu','testtime','addMultiSignalUsingExpr')
rdir <- here('hca','simu','testtime','result','addsignal', 'EM_NOT_centered')
pdir <- here('hca','simu','testtime','plot','addsignal')
Res <- readRDS(paste0(rdir, '/1.rds'))
fit <- t(sapply(names(Res$fdr[Res$fdr < 0.05]), function(g){
  tmp <- get_population_fit(Res, variable = NA, gene = g)
}))

## clustering
set.seed(12345)
clu <- kmeans(fit,10,iter.max = 1000)$cluster
clu <- sort(clu)
saveRDS(clu, paste0(rdir, '/cluster.rds'))

### use population levels
#### average the population fit for genes within a same cluster
fitClu <- sapply(1:max(clu), function(i){
  colMeans(fit[names(clu[clu==i]), 1:(ncol(fit)),drop = FALSE])
})
colnames(fitClu) <- paste0('cluster', unique(clu), '(',table(clu), ')')

##### plot population fit on original scale
pd <- reshape2::melt(fitClu)
colnames(pd) <- c('pseudotime', 'cluster', 'expr')
saveRDS(pd, paste0(pdir, '/clu_based_on_population_pattern_pd.rds'))

library(ggplot2)
cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
  
pdf(paste0(pdir, '/clu_based_on_population_pattern.pdf'), width = 4.5, height = 3.2)  
print(ggplot(data= pd, aes(x = pseudotime, y = expr, group = cluster, color = cluster)) + 
  geom_smooth(size = 1) +
  theme_classic() +
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(11, 'Set1')))(max(clu))) +
  xlab('Pseudotime') +
  ylab('Expression'))
dev.off()

##### plot population fit on standadized scale
fitClu.scale = t(scalematrix(t(fitClu)))
pd <- reshape2::melt(fitClu.scale)
colnames(pd) <- c('pseudotime', 'cluster', 'expr')
cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
saveRDS(pd, paste0(pdir, '/clu_based_on_population_pattern_scaled_pd.rds'))  
library(ggplot2)
pdf(paste0(pdir, '/clu_based_on_population_pattern_scaled.pdf'), width = 4.5, height = 3.2)  
print(ggplot(data= pd, aes(x = pseudotime, y = expr, group = cluster, color = cluster)) + 
  geom_smooth(size = 1) +
  theme_classic() +
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(11, 'Set1')))(max(clu))) +
  xlab('Pseudotime') +
  ylab('Expression'))
dev.off()

### use predict.values
pred <- Res$predict.values[names(Res$fdr[Res$fdr < 0.05]), names(Res$pseudotime)]
agg <- t(sapply(unique(clu), function(i) colMeans(pred[clu == i, ,drop=F])))
agg <- agg[, names(Res$pseudotime)]
rownames(agg) <- paste0('cluster', unique(clu), '(',table(clu), ')')

pd <- reshape2::melt(agg)
pd$pseudotime = Res$pseudotime[pd[,2]]
colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
saveRDS(pd, paste0(pdir, '/gene_cluster_pattern_pd.rds'))

pdf(paste0(pdir, '/gene_cluster_pattern.pdf'), width = 4.5, height = 3.2)  
print(ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
  geom_smooth(size = 2) +
  theme_classic() +
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(11, 'Set1')))(max(clu))) +
  xlab('Pseudotime') +
  ylab('Averaged Fitted Centered Expression'))
dev.off()



