rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
n.permute <- 10000
source("/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/function/01_function.R")
plotdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/tree_variability/2latent_clu10/plot/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/tree_variability/2latent_clu10/res/'
dir.create(rdir, recursive = T)
dir.create(plotdir, recursive = T)
# --------------------------------------------------------------
# input: seurat integrated object including:
# low dim reduction: umap, pca, or phate
# celltype: a dataframe, col 1 is cell name, col 2 is cell type (at least for the cells with origin cell type), col 3 is sample name
# origin: the origin cell type
# --------------------------------------------------------------
# read in data
pca <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/scvi/data/latent.rds')
colnames(pca) = paste0()
a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/ct/sc.rds')
ct = data.frame(cell = names(a), celltype = a, sample = sapply(names(a), function(i) sub(':.*', '', i)), stringsAsFactors = FALSE)
  
# permutation 
res = infer_tree_structure(pca = pca, ct = ct, origin.celltype = 'HSC', plotdir = plotdir, xlab='scVI latent 1', ylab = 'scVI latent 2', number.cluster=10, pcadim=NULL, original = F)
table(res$clusterid,ct[match(names(res$clusterid),ct[,1]),2])
saveRDS(res, paste0(rdir, 'infer_tree_structure_res.rds'))
png(paste0(plotdir, 'mcl.png'), width=900,height=800, res = 200)
plotmclust(res, cell_point_size = 0.1, x.lab = 'scVI latent 1', y.lab = 'scVI latent 2')
dev.off()
pdf(paste0(plotdir, 'mcl.pdf'), width=6,height=5)
plotmclust(res, cell_point_size = 0.1, x.lab = 'scVI latent 1', y.lab = 'scVI latent 2')
dev.off()
result <- evaluate_uncertainty(res, n.permute)
saveRDS(result, paste0(rdir, 'result.rds'))
for (i in names(result)){
  write.csv(result[[i]], paste0(rdir, i, '.csv'), row.names = T)
}

