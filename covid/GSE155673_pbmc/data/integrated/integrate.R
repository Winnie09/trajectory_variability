library(Matrix)
library(Seurat)
d <- readRDS('/home-4/zji4@jhu.edu/scratch/covid/data/GSE155673_pbmc/proc/count.rds')
d <- d[rowMeans(d > 0) >= 0.01,]
d <- d[!grepl('^MT-',rownames(d)),]
p <- sub(paste0(':.*'),'',colnames(d))

dl <- list()
for (i in unique(p)) {
  dl[[i]] <- CreateSeuratObject(d[,p==i])
  dl[[i]] <- NormalizeData(dl[[i]], verbose = FALSE)
  dl[[i]] <- FindVariableFeatures(dl[[i]], verbose = FALSE)
}

anchors <- FindIntegrationAnchors(object.list = dl,k.filter=NA)
integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:10,seed.use=NULL)
integrated <- FindNeighbors(integrated, reduction='umap',dims=1:2)
integrated <- FindClusters(integrated,resolution=1.2)

saveRDS(integrated,file='/home-4/zji4@jhu.edu/scratch/diffpt/covid/GSE155673_pbmc/integrate.rds')

