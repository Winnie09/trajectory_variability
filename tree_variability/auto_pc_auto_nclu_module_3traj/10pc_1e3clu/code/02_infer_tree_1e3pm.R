rm(list=ls())
library(ggplot2)
## library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
n.permute <- 1e3
max.clunum <- 1e3
nclu <- 1e3
## setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
setwd('/home/whou10/scratch16/whou10/trajectory_variability/')
source("function/01_function.R")
plotdir <- 'tree_variability/auto_pc_auto_nclu_module_3traj/10pc_1e3clu/plot/'
rdir <- 'tree_variability/auto_pc_auto_nclu_module_3traj/10pc_1e3clu/result/'

dir.create(plotdir, recursive = T)
dir.create(rdir, recursive = T)

# read in data
# umap = readRDS('hca/data/HCA/proc/integrate/ser/umap.rds')
umap = readRDS('hca/data/proc/integrate/ser/umap.rds')
pca <- as.matrix(umap@reductions$pca@cell.embeddings)
# a = readRDS('hca/data/HCA/proc/ct/sc.rds')
a = readRDS('hca/data/proc/ct/sc.rds')
ct = data.frame(cell = names(a), celltype = a, sample = sapply(names(a), function(i) sub(':.*', '', i)), stringsAsFactors = FALSE)

### infer tree structure
source('package/Lamian/R/infer_tree_structure.R')
pca = pca
cellanno = ct[,c(1,3,2)]
expression = umap@assays$RNA@counts
origin.marker = c('CD14')
origin.celltype = 'HSC'
number.cluster = nclu
plotdir = plotdir
xlab='Principal component 1'
ylab = 'Principal component 2'
max.clunum = nclu
res <- infer_tree_structure(pca = pca,
                                        cellanno = ct[,c(1,3,2)],
                                        expression = umap@assays$RNA@counts,
                                        origin.marker = c('CD14'),
                                        origin.celltype = 'HSC',
                                        number.cluster = nclu,
                                        plotdir = plotdir,
                                        xlab='Principal component 1', 
                                        ylab = 'Principal component 2',
                                        max.clunum = nclu,
                                        kmeans.seed = 7)
saveRDS(res, paste0(rdir, 'infer_tree_structure_res.rds'))
png(paste0(plotdir, 'mcl.png'), width=900,height=800, res = 200)
plotmclust(res, cell_point_size = 0.1, x.lab = 'Principal component 1', y.lab = 'Principal component 2',)
dev.off()
pdf(paste0(plotdir, 'mcl.pdf'), width=6,height=5)
plotmclust(res, cell_point_size = 0.1, x.lab = 'Principal component 1', y.lab = 'Principal component 2')
dev.off()
### evaluate branch uncertainty
result <- evaluate_uncertainty(res, n.permute)
saveRDS(result, paste0(rdir, 'result.rds'))
for (i in names(result)){
  write.csv(result[[i]], paste0(rdir, i, '.csv'), row.names = T)
}



