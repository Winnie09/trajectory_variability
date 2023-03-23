library(Seurat)
library(Matrix)
library(ggplot2)
library(umap)
d <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/ss/CD8integrate.rds')
de <- d$pca@cell.embeddings
d <- d$RNA@data

rownames(d) <- sub(':.*','',rownames(d))
#set.seed(1234)
set.seed(123)
k <- kmeans(de[,1:2],2)$cluster

library(TSCAN)
ord <- TSCANorder(exprmclust(t(de[,1:2]),cluster=k,reduce=F),orderonly=T)
saveRDS(ord,file='/home-4/zji4@jhu.edu/scratch/diffpt/covid/ss/order.rds')

sapply(c('CCR7','TCF7','NKG7','GZMB','IL7R','GZMH'),function(g) {
  cor(d[g,ord],1:length(ord))  
})

