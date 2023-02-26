pdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/slingshot/plot/'
xlab = 'Principal component 1'
ylab = 'Principal component 2'
library(Seurat)
# umap <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/integrate/ser/umap.rds')
umap = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/integrate/ser/umap.rds')
pca <- as.matrix(umap@reductions$pca@cell.embeddings)
pca <- pca[,1:2]
library(slingshot)
sce <- SingleCellExperiment(assays = List(counts = umap@assays$RNA@counts))
reducedDims(sce) <- SimpleList(PCA = pca)

### ----------
### clustering
### ----------
library(mclust, quietly = TRUE)
cl1 <- Mclust(pca)$classification
colData(sce)$GMM <- cl1

# ## plot mclust
# library(RColorBrewer)
# png(paste0(pdir, 'mclust_default.png'), width = 800, height = 800, res = 200)
# plot(pca, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1, xlab = xlab, ylab = ylab)
# dev.off()

cl2 <- kmeans(pca, centers = 4)$cluster
colData(sce)$kmeans <- cl2

# ## plot kmeans 
# png(paste0(pdir, 'kmeans_default.png'), width = 800, height = 800, res = 200)
# plot(pca, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1, xlab = xlab, ylab = ylab)
# dev.off()

### ----------------
### infer pseudotime
### ----------------
system.time({
  sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
})
summary(sce$slingPseudotime_1)
# Using full covariance matrix
#    user  system elapsed 
# 852.818   2.187 857.102 

## plot pseudotime
library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(nrow(pca))
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=nrow(pca))]

png(paste0(pdir, 'pseudotime.png'), width = 800, height = 800, res = 200)
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, xlab = xlab, ylab = ylab)
lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()

png(paste0(pdir, 'cluster.png'), width = 800, height = 800, res = 200)
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1, xlab = xlab, ylab = ylab)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
dev.off()

