library(here)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(grid)
setwd(here())
path = 'lymph'
rdir <- ddir <- paste0('hca/real/testtime/result/', path)
pdir  <- paste0('hca/real/testtime/plot/', path)
dir.create(pdir, recursive = T)
source('function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')


Res <- readRDS(paste0(ddir, '/testtime_res.rds'))

## get the differential genes
res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange, pvalue = Res$pvalue)
res <- res[order(res[,1], -res[,2]), ]

## get the population fit of the differential genes
fit <- t(sapply(names(Res$fdr[Res$fdr < 0.05]), function(g){
  tmp <- get_population_fit(Res, variable = NA, gene = g)
}))
saveRDS(fit, paste0(rdir, '/population_fit.rds')) 

fit.scale = scalematrix(fit) 

## clustering
set.seed(12345)
clu <- kmeans(fit.scale, 5, iter.max = 1000)$cluster
clu <- sort(clu)
table(clu)
saveRDS(clu, paste0(rdir, '/cluster.rds'))
  
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

  
## read in the population fit and clustering  if not in the cache
fit <- readRDS(paste0(rdir, '/population_fit.rds'))
clu <- readRDS(paste0(rdir, '/cluster.rds'))
  

res = res[res[, 1] < 0.05,]
pt = Res$pseudotime
cellanno = Res$cellanno

## standardize the values
fit.scale = scalematrix(fit)
res$clu = clu[rownames(res)]
res$cor <-
  sapply(rownames(res), function(i)
    cor(fit.scale[i, ], seq(1, ncol(fit.scale))))
fit.scale <- fit.scale[rownames(res)[order(res$clu, res$cor)],]

expr.ori <- Res$expr.ori
expr.ori <- expr.ori[rownames(fit.scale), colnames(fit.scale)]
expr.ori.scale <- scalematrix(expr.ori)

fit.scale[fit.scale > quantile(as.vector(fit.scale), 0.99)] <-
  quantile(as.vector(fit.scale), 0.99)
fit.scale[fit.scale < quantile(as.vector(fit.scale), 0.01)] <-
  quantile(as.vector(fit.scale), 0.01)
expr.ori.scale[expr.ori.scale > quantile(as.vector(expr.ori.scale), 0.99)] <-
  quantile(as.vector(expr.ori.scale), 0.99)
expr.ori.scale[expr.ori.scale < quantile(as.vector(expr.ori.scale), 0.01)] <-
  quantile(as.vector(expr.ori.scale), 0.01)
## plot
### annotate rows and columns
colann <-
  data.frame(sample = cellanno[match(cellanno[, 1], colnames(fit.scale)), 2],
             pseudotime = pt,
             stringsAsFactors = F)
rownames(colann) = colnames(fit.scale)
rowann = data.frame(cluster = as.character(clu),
                    stringsAsFactors = F)
rownames(rowann) = names(clu)

#### define colors
col.clu = colorRampPalette(brewer.pal(8, 'Set1'))(max(clu))
names(col.clu) = unique(clu)
col.sample = colorRampPalette(rev(brewer.pal(n = 8, name = "Set1")))(length(unique(colann$sample)))
names(col.sample) = unique(colann$sample)
col.pseudotime = colorRampPalette(rev(brewer.pal(n = 9, name = "YlGnBu")))(length(unique(colann$pseudotime)))
names(col.pseudotime) = unique(colann$pseudotime)

png(
  paste0(pdir, '/hm_kmeans_population_fit_scale.png'),
  width = 2100,
  height = 3200,
  res = 300
)
cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
pheatmap(
  fit.scale,
  cluster_rows = F,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = cpl,
  annotation_row = rowann,
  annotation_col = colann,
  annotation_colors = list(cluster = col.clu,
                           sample = col.sample),
  cellwidth = 250 / ncol(fit.scale),
  cellheight = 450 / nrow(fit.scale),
  border_color = NA
)

setHook("grid.newpage", NULL, "replace")
grid.text("Pseudotime", y = -0.07, gp = gpar(fontsize = 16))
grid.text(
  "Differential Genes",
  x = -0.07,
  rot = 90,
  gp = gpar(fontsize = 16)
)
dev.off()

png(
  paste0(pdir, '/hm_kmeans_original_expression_scale.png'),
  width = 2100,
  height = 3200,
  res = 300
)
cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
pheatmap(
  expr.ori.scale,
  cluster_rows = F,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = cpl,
  annotation_row = rowann,
  annotation_col = colann,
  annotation_colors = list(cluster = col.clu,
                           sample = col.sample),
  cellwidth = 250 / ncol(fit.scale),
  cellheight = 450 / nrow(fit.scale),
  border_color = NA
)
dev.off()
}

