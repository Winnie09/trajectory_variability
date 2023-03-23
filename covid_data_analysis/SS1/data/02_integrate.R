library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)

npcs = 20
nfeatures = 2000
set.seed(12345)
u <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/ss/count.rds')
sample <- sub(paste0(':.*'),'',colnames(u))
u <- u[rowMeans(u > 0) >= 0.01,]
u <- u[!grepl('^MT-',rownames(u)),]

## 1. harmony
u <- CreateSeuratObject(counts = u, project = "covid")
u <- NormalizeData(u)
u <- FindVariableFeatures(u,selection.method = "vst", nfeatures = nfeatures)
u <- ScaleData(u,verbose = FALSE)
u <- RunPCA(u,npcs = npcs, verbose = FALSE)

u@meta.data$sample <- sample
u <- RunHarmony(u, 'sample')
u <- RunUMAP(u, reduction = "harmony", dims = 1:npcs)
u <- FindNeighbors(u,reduction = "harmony", dims = 1:npcs)
u <- FindClusters(u)
saveRDS(u,'/home-4/zji4@jhu.edu/scratch/diffpt/covid/ss/integrate.rds')


