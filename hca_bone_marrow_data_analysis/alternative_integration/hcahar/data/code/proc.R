library(Seurat)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/integrate/harmony.rds')
pca <- d$harmony@cell.embeddings
u <- d$umap@cell.embeddings
clu <- Idents(d)
n <- names(clu)
clu <- as.numeric(clu) + 1
names(clu) <- n
saveRDS(pca,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/proc/pca.rds')
saveRDS(u,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/proc/umap.rds')
saveRDS(clu,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/proc/cluster.rds')

