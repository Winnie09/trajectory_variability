rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
n.permute <- 10000
max.clunum <- 50
source("/Users/wenpinhou/Dropbox/trajectory_variability/function/01_function.R")
plotdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/auto_pc_auto_nclu_module_3traj/plot/'
rdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/auto_pc_auto_nclu_module_3traj/result/'
# --------------------------------------------------------------
# input: seurat integrated object including:
# low dim reduction: umap, pca, or phate
# celltype: a dataframe, col 1 is cell name, col 2 is cell type (at least for the cells with origin cell type), col 3 is sample name
# origin: the origin cell type
# --------------------------------------------------------------
# read in data
umap = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/integrate/ser/umap.rds')
pca <- as.matrix(umap@reductions$pca@cell.embeddings)
a = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/ct/sc.rds')
ct = data.frame(cell = names(a), celltype = a, sample = sapply(names(a), function(i) sub(':.*', '', i)), stringsAsFactors = FALSE)
  
# permutation 
res = infer_tree_structure(pca = pca, ct = ct, origin.celltype = 'HSC', plotdir = plotdir, xlab='Principal component 1', ylab = 'Principal component 2', original = T)
saveRDS(res, paste0(rdir, 'infer_tree_structure_res.rds'))
png(paste0(plotdir, 'mcl.png'), width=900,height=800, res = 200)
plotmclust(res, cell_point_size = 0.1, x.lab = 'Principal component 1', y.lab = 'Principal component 2')
dev.off()
pdf(paste0(plotdir, 'mcl.pdf'), width=6,height=5)
plotmclust(res, cell_point_size = 0.1, x.lab = 'Principal component 1', y.lab = 'Principal component 2')
dev.off()
result <- evaluate_uncertainty(res, n.permute, design=data.frame(group=ifelse(paste0('BM',1:8) %in% c('BM1', 'BM2', 'BM5', 'BM6'),'G1','G2'),row.names = paste0('BM',1:8)))
saveRDS(result, paste0(rdir, 'result.rds'))

for (i in names(result)){
  write.csv(result[[i]], paste0(rdir, i, '.csv'), row.names = T)
}

## subsample cells, and then redo infer tree structure
# ---------------
# for all samples
# ---------------
for (rm.perc in seq(0.1, 0.8, 0.1)){
  print(rm.perc)
  plotdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmall/', rm.perc, '/plot/')
  rdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmall/', rm.perc, '/result/')
  dir.create(plotdir, recursive = T)
  dir.create(rdir, recursive = T)
  
  selectcell = res$order[[2]]  ## get branch 5,1 cells
  set.seed(12345)
  rmcell = sample(selectcell, rm.perc*length(selectcell)) 
  subset.cell = setdiff(rownames(pca), rmcell) ## remove a percentage of cells from branch 5,1

  pdf(paste0(plotdir, 'mcl.pdf'), width=6,height=5)
  print(plotmclust(res, cell_point_size = 0.1, x.lab = 'Principal component 1', y.lab = 'Principal component 2', subset.cell = subset.cell))
  dev.off()
  result <- evaluate_uncertainty(res, n.permute, subset.cell = subset.cell)
  saveRDS(result, paste0(rdir, 'result.rds'))
  for (i in names(result)){
    write.csv(result[[i]], paste0(rdir, i, '.csv'), row.names = T)
  }
}
  
# ----------------------------
# for some samples: BM1,2,5,6
# ----------------------------
for (rm.perc in seq(0.1, 0.8, 0.1)){
  print(rm.perc)
  plotdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256/', rm.perc, '/plot/')
  rdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256/', rm.perc, '/result/')
  dir.create(plotdir, recursive = T)
  dir.create(rdir, recursive = T)
  
  selectcell = res$order[[2]]
  selectcell = selectcell[ct[selectcell, 'sample'] %in% c('BM1', 'BM2', 'BM5', 'BM6')]
  set.seed(12345)
  rmcell = sample(selectcell, rm.perc*length(selectcell))
  subset.cell = setdiff(rownames(pca), rmcell)

  # pdf(paste0(plotdir, 'mcl.pdf'), width=6,height=5)
  # print(plotmclust(res, cell_point_size = 0.1, x.lab = 'PC1', y.lab = 'PC2', subset.cell = subset.cell))
  # dev.off()
  result <- evaluate_uncertainty(res, n.permute, subset.cell = subset.cell,design=data.frame(group=ifelse(paste0('BM',1:8) %in% c('BM1', 'BM2', 'BM5', 'BM6'),'G1','G2'),row.names = paste0('BM',1:8)))
  saveRDS(result, paste0(rdir, 'result.rds'))
  for (i in names(result)){
    write.csv(result[[i]], paste0(rdir, i, '.csv'), row.names = T)
  }
}
  
