setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/proc')
library(Seurat)
u = readRDS('./integrate/umap.rds')
u = u@reductions$umap@cell.embeddings  
saveRDS(u,'./integrate/umap/umap.rds')
sample = gsub(':.*','',rownames(u))
cluster = readRDS('./cluster/cluster.rds')
ct = readRDS('./ct/sc.rds')

pd = data.frame(cell = rownames(u),umap1 = u[,1], umap2 = u[,2], sample = sample, cluster = as.factor(cluster[rownames(u)]), ct = ct[rownames(u)])

library(ggplot2)
p <- ggplot() + geom_point(data = pd, aes(x=umap1, y=umap2, color=sample),alpha=0.5,size=0.2) + theme_classic() + xlab('UMAP1')+ylab('UMAP2')
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/umap_sample.png',
        p,
        width=5.5,height=4,dpi=300)


p <- ggplot() + geom_point(data = pd, aes(x=umap1, y=umap2, color=ct),alpha=0.5,size=0.2) + theme_classic() + xlab('UMAP1')+ylab('UMAP2') +
  guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/umap_ct.png',
        p,
        width=5.5,height=4,dpi=300)

p <- ggplot() + geom_point(data = pd, aes(x=umap1, y=umap2, color=ct),alpha=0.5,size=0.1) + theme_classic() + xlab('UMAP1')+ylab('UMAP2') +
  guides(color = guide_legend(override.aes = list(size = 3,alpha=1))) + facet_wrap(~ct)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/umap_ct_facet.png',
        p,
        width=8,height=6,dpi=300)


cm = aggregate(u,list(cluster),mean)
colnames(cm) <- c('cluster','x','y')
p <- ggplot() + geom_point(data=pd, aes(x=umap1,y=umap2,col=cluster),alpha=0.5,size=0.2) + 
  geom_text(data=cm,aes(x=x,y=y,label=cluster),size=3) + theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'none')
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/umap_cluster.png',
        p,
        width=5.5,height=4,dpi=300)
