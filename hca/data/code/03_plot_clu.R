setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/')
library(Seurat)
set.seed(12345)
integrated <- readRDS('./proc/integrate/umap.rds')
resolution = 0.1
d <- FindNeighbors(integrated, reduction='umap',dims=1:2)
d <- FindClusters(d,resolution=resolution)
clu = Idents(d)
str(clu)
u = readRDS('./proc/integrate/umap/umap.rds')
# plot(u,col=clu)
saveRDS(clu, './proc/cluster/resolution0.1.rds')

library(TSCAN)
set.seed(12345)
em <- exprmclust(t(u),reduce=F,cluster = clu)
# plotmclust(em)

ct = readRDS('./proc/ct/sc.rds')
ct = ct[rownames(u)]
clucenter <- aggregate(u,list(clu),mean)
colnames(clucenter) <- c('cluster','x','y')

library(igraph)
edgelist <- matrix(as.numeric(as_edgelist(em$MSTtree)),ncol=2)
ap <- data.frame(x=clucenter[edgelist[,1],'x'],y=clucenter[edgelist[,1],'y'],xend=clucenter[edgelist[,2],'x'],yend=clucenter[edgelist[,2],'y'])

pd1 = data.frame(x=u[,1],y=u[,2],ct=ct,clu = clu)
library(ggplot2)
library(gridExtra)
p1<- ggplot() + geom_point(data=pd1,aes(x=x,y=y,col=ct),alpha=0.5,size=0.2) + 
  # geom_segment(data=ap,aes(x=x,y=y,xend=xend,yend=yend),size=1,color='brown') + 
  geom_text(data=clucenter,aes(x=x,y=y,label=cluster),size=7) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme_classic() + theme(legend.title = element_blank()) +
  xlab('UMAP1') + ylab('UMAP2')
p2<- ggplot() + geom_point(data=pd1,aes(x=x,y=y,col=clu),alpha=0.5,size=0.2) + 
  geom_text(data=clucenter,aes(x=x,y=y,label=cluster),size=7) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme_classic() + theme(legend.title = element_blank()) +
  xlab('UMAP1') + ylab('UMAP2')
ggsave('./plot/pseudotime.png',
       grid.arrange(p1,p2,nrow=1),
       width=12,height=5,
       dpi=200)

