library(Seurat)
library(TSCAN)
d <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/GSE155673_pbmc/integrate.rds')
gl <- sapply(c('CD3D:','CD4:','CD8A:','TCF7:'),function(i) grep(i,rownames(d),value=T))
expr <- AverageExpression(d,features = gl)$RNA
den1 <- density(unlist(expr[1,]))
den3 <- density(unlist(expr[3,]))
clu <- colnames(expr[,expr[1,] > den1$x[sum(rle(diff(den1$y) > 0)$lengths[1:2])] & expr[3,] > den3$x[sum(rle(diff(den3$y) > 0)$lengths[1:2])]])
k <- Idents(d)
n <- names(k)
k <- as.character(k)
names(k) <- n
cell <- names(k)[which(k %in% clu)]
saveRDS(cell,file='/home-4/zji4@jhu.edu/scratch/diffpt/covid/GSE155673_pbmc/CD8cell.rds')

