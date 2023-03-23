library(Seurat)
library(Matrix)
library(pheatmap)
d <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/ss/integrate.rds')
e <- d$RNA@data
rownames(e) <- sub(':.*','',rownames(e))
clu <- Idents(d)
m <- aggregate(as.matrix(t(e[c('CD3D','CD4','CD8A','CD14','CD19','CCR7','TCF7','NKG7','GZMB'),])),list(clu),mean)
rownames(m) <- m[,1]
m <- as.matrix(m[,-1])
#pheatmap(m,scales='row')
# pdf('clu.pdf')
# DimPlot(d,label = T)
# dev.off()
den1 <- density(unlist(m[,'CD3D']),bw = 'SJ')
den3 <- density(unlist(m[,'CD8A']),bw = 'SJ')

clu <- rownames(m)[m[,'CD3D'] > 1 & m[,'CD8A'] > 0.5]

print(clu)
k <- Idents(d)
n <- names(k)
k <- as.character(k)
names(k) <- n
cell <- names(k)[which(k %in% clu)]
tab <- table(sub(':.*','',cell))
n <- names(tab)[tab > 50]
cell <- cell[sub(':.*','',cell) %in% n]
saveRDS(cell,file='/home-4/zji4@jhu.edu/scratch/diffpt/covid/ss/CD8cell.rds')


