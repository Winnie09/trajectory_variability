d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/data/latent.rds')
sc <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/ct/sc_cluster.rds')
clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/cluster/cluster.rds')
ct <- sc[clu[rownames(d)]]
names(ct) <- names(clu)

lct <- ct
lct[ct %in% c('Bcell','CD4Tcell','CD8Tcell','NKcell')] <- 'immune'
lct[ct %in% c('Ery','MEP')] <- 'ery'
lct[ct %in% c('HSC','LMPP','GMP')] <- 'hsc'

set.seed(1)
clu <- kmeans(d,7)$cluster
table(lct,clu)
library(TSCAN)
mcl <- exprmclust(t(d),cluster=clu)
#plot(mcl$MSTtree)
library(igraph)
library(RColorBrewer)

mcl$MSTtree <- mcl$MSTtree - edges('1','2') - edges('1','4') - edges('1','2') + edges('2','7') + edges('7','1') + edges('1','5') + edges('2','6') + edges('6','4') + edges('2','3')

saveRDS(mcl, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/data/mcl.rds')
#mcl$MSTtree <- minimum.spanning.tree(graph.data.frame(cbind(3,c(1,2,4)),directed = F))
o <- TSCANorder(mcl,orderonly = T,listbranch = T,startcluster = 2)
sapply(o,function(i) table(ct[i]),simplify = F)
saveRDS(o,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/data/order_originalNames.rds')

#2-7-1-5, 2-3, 2-6-4
names(o) <- c('lymphoid','myeloid','erythroid')
saveRDS(o,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/data/order.rds')

med <- aggregate(d,list(clu),median)
library(ggplot2)
cn1 <- 'V1'
cn2 <- 'V2'

ldm <- as_edgelist(mcl$MSTtree)
ld <- do.call(rbind,sapply(1:nrow(ldm),function(i) {
  i=as.numeric(ldm[i,])
  data.frame(x=med[i[1],cn1],xend=med[i[2],cn1],y=med[i[1],cn2],yend=med[i[2],cn2])
},simplify = F))
library(scattermore)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/plot/d1d2.pdf',width=3.1,height=3.2)
ggplot() + geom_scattermore(data=data.frame(x=d[,1],y=d[,2],ct=as.character(clu)),aes(x=x,y=y,col=ct),size=0.1, stroke = 0) + geom_text(data=med,aes(x=V1,y=V2,label=Group.1),size=5) + geom_segment(data=ld,aes(x=x,xend=xend,y=y,yend=yend))  + theme(legend.position = 'none', axis.text = element_text(size = 16)) + xlab('scVI latent 1') + ylab('scVI latent 2') + scale_color_brewer(palette = 'Dark2') 
dev.off()


cn1 <- 'V1'
cn2 <- 'V3'
ldm <- as_edgelist(mcl$MSTtree)
ld <- do.call(rbind,sapply(1:nrow(ldm),function(i) {
  i=as.numeric(ldm[i,])
  data.frame(x=med[i[1],cn1],xend=med[i[2],cn1],y=med[i[1],cn2],yend=med[i[2],cn2])
},simplify = F))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/plot/d1d3.pdf',width=4,height=4)
ggplot() + geom_point(data=data.frame(x=d[,1],y=d[,3],ct=as.character(clu)),aes(x=x,y=y,col=ct),size=0.1) + geom_text(data=med,aes(x=V1,y=V3,label=Group.1),size=5) + geom_segment(data=ld,aes(x=x,xend=xend,y=y,yend=yend)) + theme_classic() + theme(legend.position = 'none') + xlab('Dim1') + ylab('Dim3')
dev.off()


cn1 <- 'V2'
cn2 <- 'V3'
ldm <- as_edgelist(mcl$MSTtree)
ld <- do.call(rbind,sapply(1:nrow(ldm),function(i) {
  i=as.numeric(ldm[i,])
  data.frame(x=med[i[1],cn1],xend=med[i[2],cn1],y=med[i[1],cn2],yend=med[i[2],cn2])
},simplify = F))

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/plot/d2d3.pdf',width=4,height=4)
ggplot() + geom_point(data=data.frame(x=d[,2],y=d[,3],ct=as.character(clu)),aes(x=x,y=y,col=ct),size=0.1) + geom_text(data=med,aes(x=V2,y=V3,label=Group.1),size=5) + geom_segment(data=ld,aes(x=x,xend=xend,y=y,yend=yend)) + theme_classic() + theme(legend.position = 'none') + xlab('Dim2') + ylab('Dim3')
dev.off()




