rm(list = ls())
library(ggplot2)
library(reshape2)
library(RColorBrewer)
suppressMessages(library(igraph))
library(parallel)
library(splines)
library(viridis)
library(pheatmap)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

dataType = 4
ddir <- here('hca','data','simu','testtime','addMultiSignalUsingExpr')
rdir <- here('hca','simu','testtime','old_addMultiSignalUsingExpr_use_all_othgene_10cluster','result','addsignal', 'cluster', dataType)
pdir <- here('hca','simu','testtime','old_addMultiSignalUsingExpr_use_all_othgene_10cluster','plot','addsignal', dataType)
Res <- readRDS(here('hca','simu','testtime','old_addMultiSignalUsingExpr_use_all_othgene_10cluster','result','addsignal', 'EM_NOT_centered',paste0(dataType,'.rds')))

clu <- readRDS(paste0(rdir, '/cluster.rds'))
fit <- readRDS(paste0(rdir, '/population_fit.rds'))

res = data.frame(
fdr = Res$fdr,
foldchange = Res$foldchange,
pvalue = Res$pvalue
)
res = res[res[, 1] < 0.05, ]
res = res[order(res[, 1], res[, 2]), ]
str(res)
pt = Res$pseudotime
cellanno = Res$cellanno

## standardize the values
fit.scale = scalematrix(fit)
res$clu = clu[rownames(res)]
res$cor <- sapply(rownames(res), function(i) cor(fit.scale[i,], seq(1,ncol(fit.scale))))
fit.scale <- fit.scale[rownames(res)[order(res$clu, res$cor)], ]
# write.csv(res, paste0(rdir, '/differential_genes_clu.csv'))

## -------
## plot
## -------
# annotate rows and columns
colann <-
data.frame(sample = cellanno[match(cellanno[, 1], colnames(fit.scale)), 2],
           stringsAsFactors = F)
rownames(colann) = colnames(fit.scale)
rowann = data.frame(cluster = as.character(clu), stringsAsFactors = F)
rownames(rowann) = names(clu)

# define colors
library(pheatmap)
col.clu = colorRampPalette(brewer.pal(8, 'Set1'))(max(clu))
names(col.clu) = unique(clu)
col.sample = colorRampPalette(rev(brewer.pal(n = 8, name = "Set1")))(length(unique(colann$sample)))
names(col.sample) = unique(colann$sample)

png(
  paste0(pdir, '/hm_kmeans_scale.png'),
  width = 2100,
  height = 3200,
  res = 300
)
cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
# cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MAC1
pheatmap(
fit.scale,
cluster_rows = F,
cluster_cols = FALSE,
show_rownames = FALSE,
show_colnames = FALSE,
color = cpl,
annotation_row = rowann,
annotation_col = colann,
annotation_colors = list(
  cluster = col.clu,
  sample = col.sample
),
cellwidth = 241.24 / ncol(fit.scale),
# cellheight = 467.35 / nrow(fit.scale),
cellheight = 450 / nrow(fit.scale),
border_color = NA
)
dev.off()


