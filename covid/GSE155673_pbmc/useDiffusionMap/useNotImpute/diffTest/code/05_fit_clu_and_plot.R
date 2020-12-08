library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))

type = 'tex_testtime'
rdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type)
pdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','plot', type)
dir.create(pdir)
Res <- readRDS(here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type, 'res.rds'))
fit <- t(sapply(names(Res$fdr[Res$fdr < 0.05]), function(g){
  tmp <- get_population_fit(Res, variable = NA, gene = g)
}))
saveRDS(fit, paste0(rdir, '/population_fit.rds')) ##########

fit.scale = scalematrix(fit) #############
## -----------
## clustering
## -----------
set.seed(12345)
clu <- kmeans(fit.scale, 5, iter.max = 1000)$cluster
clu <- sort(clu)
table(clu)
saveRDS(clu, paste0(rdir, '/cluster.rds'))

### use population levels
#### average the population fit for genes within a same cluster
fitClu <- sapply(1:max(clu), function(i){
  colMeans(fit[names(clu[clu==i]), 1:(ncol(fit)),drop = FALSE])
})
colnames(fitClu) <- paste0('cluster', unique(clu), '(',table(clu), ')')

##### plot population fit on original scale
pd <- reshape2::melt(fitClu)
colnames(pd) <- c('cell', 'cluster', 'expr')
pd$pseudotime <- Res$pseudotime[pd[,1]]
cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
saveRDS(pd, paste0(pdir, '/clu_based_on_population_pattern_pd.rds'))
pdf(paste0(pdir, '/clu_based_on_population_pattern.pdf'), width = 4.5, height = 3)  
print(ggplot(data= pd, aes(x = pseudotime, y = expr, group = cluster, color = cluster)) + 
  geom_smooth(size = 1) +
  theme_classic() +
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, 'Set1')))(max(clu))) +
  xlab('Pseudotime') +
  ylab('Expression') +
    theme(axis.text.x = element_blank()))
dev.off()

##### plot population fit on standadized scale
fitClu.scale = t(scalematrix(t(fitClu)))
pd <- reshape2::melt(fitClu.scale)
colnames(pd) <- c('cell', 'cluster', 'expr')
pd$pseudotime <- Res$pseudotime[pd[,1]]
cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
saveRDS(pd, paste0(pdir, '/clu_based_on_population_pattern_scaled_pd.rds'))  
pdf(paste0(pdir, '/clu_based_on_population_pattern_scaled.pdf'), width = 4.5, height = 3.2)  
print(ggplot(data= pd, aes(x = pseudotime, y = expr, group = cluster, color = cluster)) + 
  geom_smooth(size = 1) +
  theme_classic() +
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, 'Set1')))(max(clu))) +
  xlab('Pseudotime') +
  ylab('Expression')) +
  theme(axis.text.x = element_blank())
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
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, 'Set1')))(max(clu))) +
  xlab('Pseudotime') +
  ylab('Averaged Fitted Expression'))
dev.off()

#### scaled 
agg.scale <- scalematrix(agg)
pd <- reshape2::melt(agg.scale)
pd$pseudotime = Res$pseudotime[pd[,2]]
colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
saveRDS(pd, paste0(pdir, '/gene_cluster_pattern_scaled_pd.rds'))

pdf(paste0(pdir, '/gene_cluster_pattern_scaled.pdf'), width = 4.5, height = 3.2)  
print(ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
  geom_smooth(size = 2) +
  theme_classic() +
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, 'Set1')))(max(clu))) +
  xlab('Pseudotime') +
  ylab('Scaled Averaged Fitted Expression'))
dev.off()

