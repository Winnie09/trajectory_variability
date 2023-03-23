library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)

npcs = 20
nfeatures = 2000

set.seed(12345)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/proc/matrix/count.rds')
p <- sub(paste0(':.*'),'',colnames(d))
d <- d[rowMeans(d > 0) > 0.01,]

## 1. harmony
d <- CreateSeuratObject(counts = d, project = "covid")
d <- NormalizeData(d)
d <- FindVariableFeatures(d,selection.method = "vst", nfeatures = nfeatures)
d <- ScaleData(d,verbose = FALSE)
d <- RunPCA(d,npcs = npcs, verbose = FALSE)

d@meta.data$p <- p
u <- RunHarmony(d, 'p')
u <- RunUMAP(u, reduction = "harmony", dims = 1:npcs)
u <- FindNeighbors(u,reduction = "harmony", dims = 1:npcs)
u <- FindClusters(u)
saveRDS(u,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/integrate/harmony.rds')


