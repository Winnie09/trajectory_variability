setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/HCA/')
set.seed(12345)
u = readRDS('./proc/integrate/umap/umap.rds')
clu = readRDS('./proc/cluster/resolution0.1.rds')
n <- names(clu)
clu = as.numeric(as.character(clu))
names(clu) = n
clu = as.factor(clu)
order = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/result/ery/order.rds')
Patient = gsub('_.*','',order$Cell)
g1 = unique(Patient[grepl('female',Patient)])
g2 = setdiff(unique(Patient), g1)
ct = readRDS('./proc/ct/sc.rds')
ct = ct[rownames(u)]
clucenter <- aggregate(u,list(clu),mean)
colnames(clucenter) <- c('cluster','x','y')
clucenter = clucenter[clucenter[,1]%in% c(7,14,8,3,0,4),]

order = data.frame(order, Patient = Patient, UMAP1 = u[order$Cell,1], UMAP2 = u[order$Cell,2], clu = clu[order$Cell], ct = ct[order$Cell], Age = sapply(Patient, function(i) strsplit(i,':')[[1]][2]), Gender = sapply(Patient, function(i) strsplit(i,':')[[1]][3]))

order$Gender = factor(as.character(order$Gender), levels = unique(order$Gender))
order$Age = factor(as.character(order$Age), levels = as.character(sort((unique(order$Age)))))
order$clu = factor(as.character(order$clu), levels = as.character(unique(order$clu)))
order$Patient = factor(as.character(order$Patient), levels = c(g1,g2))

library(ggplot2)
library(gridExtra)


p1<- ggplot() + geom_point(data=order,aes(x=UMAP1,y=UMAP2,col=ct),alpha=0.5,size=0.2) + 
  geom_text(data=clucenter,aes(x=x,y=y,label=cluster),size=7) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme_classic() + theme(legend.title = element_blank()) +
  facet_wrap(~Patient,nrow=2)+
  xlab('UMAP1') + ylab('UMAP2')
p2<- ggplot() + geom_point(data=order,aes(x=UMAP1,y=UMAP2,col=clu),alpha=0.5,size=0.2) + 
  geom_text(data=clucenter,aes(x=x,y=y,label=cluster),size=7) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme_classic() + theme(legend.title = element_blank()) +
  facet_wrap(~Patient,nrow=2)+
  xlab('UMAP1') + ylab('UMAP2')
p3<- ggplot() + geom_point(data=order,aes(x=UMAP1,y=UMAP2,col=order$Pseudotime),alpha=0.5,size=0.2) + 
  theme_classic() + theme(legend.title = element_blank()) +
  facet_wrap(~Patient,nrow=2)+
  xlab('UMAP1') + ylab('UMAP2')+
  scale_color_gradient2(low='darkblue', mid='green', high = 'yellow',midpoint=7500)

ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/plot/plot/ery_ct_clu_pseudotime.png',
       grid.arrange(p1,p2,p3, nrow=1),
       width=20,height=3,
       dpi=200)
