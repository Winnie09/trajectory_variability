rm(list = ls())
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
set.seed(12345)
pca <- readRDS('hca/data/HCA/proc/integrate/pca/pca.rds')
order = readRDS('hca/result/ery/order.rds')
Patient = gsub('_.*','',order$Cell)
g1 = unique(Patient[grepl('female',Patient)])
g2 = setdiff(unique(Patient), g1)
ct = readRDS('hca/data/HCA/proc/ct/sc.rds')
ct = ct[rownames(pca)]
order = data.frame(order, Patient = Patient, PC1 = pca[order$Cell,1], PC2 = pca[order$Cell,2])

order$Patient = factor(as.character(order$Patient), levels = c(g1,g2))

library(ggplot2)
library(gridExtra)
library(RColorBrewer)

mycolor <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(1000)
p<- ggplot() + geom_point(data=order,aes(x=PC1,y=PC2,col=order$Pseudotime),alpha=0.5,size=0.2) + 
  theme_classic() + theme(legend.title = element_blank()) +
  facet_wrap(~Patient,nrow=2)+
  xlab('Principal component 1') + ylab('Principal component 2')+
    scale_color_gradientn(colors = mycolor)
ggsave('hca/plot/plot/ery_ct_clu_pseudotime_pca.png',
       p,
       width=5.8,height=2.7,
       dpi=300)


