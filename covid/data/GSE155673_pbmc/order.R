library(Seurat)
library(TSCAN)
d <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/GSE155673_pbmc/CD8integrate.rds')
gl <- c('CCR7','TCF7','GZMK','XCL1','CX3CR1','HAVCR2')
gl <- paste0(gl,':')
gl <- sapply(gl,function(i) grep(i,rownames(d),value=T))
pdf('/home-4/zji4@jhu.edu/scratch/diffpt/covid/GSE155673_pbmc/marker.pdf')
FeaturePlot(d,gl,cols=rainbow(15)[c(12:1,15,15)])
dev.off()
u <- d$umap@cell.embeddings
cell <- which(u[,2] < -3 & u[,1] > 5)
u <- u[setdiff(1:nrow(u),cell),]
set.seed(12345)
clu <- kmeans(u,7)$cluster
pdf('/home-4/zji4@jhu.edu/scratch/diffpt/covid/GSE155673_pbmc/umap_cluster.pdf')
plot(u,col=clu)
dev.off()
#mc <- exprmclust(t(u),cluster=clu,reduce=F)
#plotmclust(mc)


p <- which(clu %in% c(2,6,3,5))
cl <- clu[p]
n <- names(cl)
cl <- as.numeric(as.factor(cl))
names(cl) <- n
mc <- exprmclust(t(u[p,]),cluster=cl,reduce=F)
ord1 <- TSCANorder(mc,orderonly=T)

p <- which(clu %in% c(2,6,7,4,1))
cl <- clu[p]
n <- names(cl)
cl <- as.numeric(as.factor(cl))
names(cl) <- n
mc <- exprmclust(t(u[p,]),cluster=cl,reduce=F)
ord2 <- TSCANorder(mc,orderonly=T)
ord <- list(temra=ord1,tex=ord2)
saveRDS(ord,file='/home-4/zji4@jhu.edu/scratch/diffpt/covid/GSE155673_pbmc/order.rds')
