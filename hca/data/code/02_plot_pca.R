setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/proc')
library(Seurat)
u = readRDS('./integrate/umap.rds')
u = u@reductions$pca@cell.embeddings  
saveRDS(u,'./integrate/pca_only.rds')
sample = gsub(':.*','',rownames(u))
cluster = readRDS('./cluster/cluster.rds')
ct = readRDS('./ct/sc.rds')
mat = readRDS('./matrix/normcount.rds')
rownames(mat) = gsub(':.*','',rownames(mat))
mat = mat[,rownames(u)]
pd = data.frame(cell = rownames(u),umap1 = u[,1], umap2 = u[,2], umap3 = u[,3], sample = sample, cluster = as.factor(cluster[rownames(u)]), ct = ct[rownames(u)])

library(ggplot2)
p <- ggplot() + geom_point(data = pd, aes(x=umap1, y=umap2, color=sample),alpha=0.5,size=0.2) + theme_classic() + xlab('PC1')+ylab('PC2')
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/pca_sample.png',
        p,
        width=5.5,height=4,dpi=300)

colv = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/colors3.rds')
colv = colv[1:length(unique(ct))]
names(colv) = (unique(ct))
p <- ggplot() + geom_point(data = pd, aes(x=umap1, y=umap2, color=ct),alpha=0.5,size=0.2) + theme_classic() + xlab('PC1')+ylab('PC2') + guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))+scale_color_manual(values=colv)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/pca_ct.png',
        p,
        width=5.5,height=4,dpi=300)
p <- ggplot() + geom_point(data = pd, aes(x=umap1, y=umap3, color=ct),alpha=0.5,size=0.2) + theme_classic() + xlab('PC1')+ylab('PC2') + guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))+scale_color_manual(values=colv)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/pca_ct13.png',
        p,
        width=5.5,height=4,dpi=300)
p <- ggplot() + geom_point(data = pd, aes(x=umap2, y=umap3, color=ct),alpha=0.5,size=0.2) + theme_classic() + xlab('PC1')+ylab('PC2') + guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))+scale_color_manual(values=colv)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/pca_ct23.png',
        p,
        width=5.5,height=4,dpi=300)
png('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/pca_ct2.png',width=700,height=700)
scatterplot3d(pd[,2:4], pch = 20, 
          color = colv[match( pd$ct,names(colv))],
          xlab='PC1',ylab='PC2',zlab='PC3')
legend("right", legend = names(colv),
      col = colv, pch = 16,cex=1)
dev.off()

cm = aggregate(u,list(cluster),mean)
colnames(cm) <- c('cluster','x','y')
p <- ggplot() + geom_point(data=pd, aes(x=umap1,y=umap2,col=cluster),alpha=0.5,size=0.2) + 
  geom_text(data=cm,aes(x=x,y=y,label=cluster),size=3) + theme_classic() + xlab('PC1')+ylab('PC2')+
  theme(legend.position = 'none')
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/plot/pca_cluster.png',
        p,
        width=5.5,height=4,dpi=300)

