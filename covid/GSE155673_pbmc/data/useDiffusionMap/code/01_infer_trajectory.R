library(destiny)
library(here)
here()
plotdir <- here('covid','data','GSE155673_pbmc','useDiffusionMap','plot')
rdir <- here('covid','data','useDiffusionMap','result')
# read in data, filtering
seu <- readRDS(here('covid', 'data','GSE155673_pbmc','integrated','CD8integrate.rds'))
expr <- as.matrix(seu@assays$RNA@data) ## librarysize-normalized  log2-transform
str(expr)
expr <- expr[rowMeans(expr>0.1)>0.01, ] ## nothing changed
str(expr) ##

# generate diffusion map
dm <- DiffusionMap(t(expr), )
str(dm)
saveRDS(dm, here('covid', 'data', 'GSE155673_pbmc','useDiffusionMap', 'result','dm.obj.rds'))
dmap <- dm@eigenvectors
rownames(dmap) <- colnames(expr)
str(dmap)
saveRDS(dmap, here('covid','data','GSE155673_pbmc','useDiffusionMap','result','dm.rds'))

# plot diffusion map
pdf(paste0(plotdir, '/dm_12dc.pdf'), width = 6, height = 5)
plot(dm,1:2,
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dm,2:3,
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dm,c(1,3),
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_3dc.pdf'), width = 6, height = 6)
plot(dm,1:3,
pch = 20) 
dev.off()

# clustering
library(RColorBrewer)
library(scattermore)
library(ggplot2)
source(here('function','01_function.R'))
number.cluster <- NA
set.seed(12345)
clu <- mykmeans(dmap[,c(1,3)], maxclunum = 50, number.cluster = number.cluster)$cluster
table(clu)
saveRDS(clu, paste0(rdir, '/cluster.rds'))
pdf(paste0(plotdir, '/cluster.pdf'), width = (0.8*max(clu))*2, height = (0.5*max(clu)))
p1 <- ggplot(data = data.frame(x = dmap[,1], y = dmap[,2], clu = as.factor(clu[rownames(dmap)])), 
             aes(x = x, y = y, color = clu)) +
  geom_scattermore()+
  scale_color_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(max(clu)))+
  theme_classic() + xlab('DM1') + ylab('DM2')
p2 <- ggplot(data = data.frame(x = dmap[,1], y = dmap[,3], clu = as.factor(clu[rownames(dmap)])), 
             aes(x = x, y = y, color = clu)) +
  geom_scattermore()+
  scale_color_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(max(clu)))+
  theme_classic() + xlab('DM1') + ylab('DM3')
gridExtra::grid.arrange(p1,p2, nrow=1)
dev.off()

# check marker genes to roughly identify cell types for pseudotime directions
# Conclusion: 1:temra 3:naive 4:Tex

allg <- c('CCR7', 'HAVCR2', 'CCL5','GZMH','GZMK', 'CX3CR1', 'TCF7', 'XCL1')
for (g in allg){
  pdf(paste0(plotdir, '/marker_genes_', g, '.pdf'), width = 3.2, height = 2.2)
  print(ggplot(data = data.frame(x = dmap[,1], y = dmap[,3], log2norm = expr[grep(paste0(g,':'), rownames(expr))[1],]), 
               aes(x = x, y = y, color = log2norm)) +
    geom_scattermore(size = 1)+
    scale_color_gradient(low = 'grey', high = 'red')+
    theme_classic() + xlab('DM1') + ylab('DM3')  +
    ggtitle(g))
  dev.off()
}

allg <- c('CCR7','HAVCR2','CCL5','GZMH')
sapply(allg,function(g) {
tapply(expr[grep(paste0(g,':'), rownames(expr))[1],],list(clu),mean)  
})
#         CCR7     HAVCR2     CCL5       GZMH
# 1 0.02145070 0.01271862 3.567711 2.44012085
# 2 0.09810904 0.06730663 3.176577 1.89181104
# 3 0.75350242 0.02624854 0.520584 0.04439525
# 4 0.04344682 0.18416064 3.219459 2.45375778

fc <- rowMeans(expr[,names(clu)[clu==1]]) - rowMeans(expr[,names(clu)[clu!=1]])

# infer pseudotime  
library(TSCAN)
mc <- exprmclust(t(dmap[,c(1,3)]),cluster=clu,reduce = T)
plot(mc$MSTtree)

ord1 <- TSCANorder(mc, orderonly = T,MSTorder = c(3,2,4))
saveRDS(ord1,here('covid','data','GSE155673_pbmc','useDiffusionMap','tex','pseudotime_order.rds'))

ord2 <- TSCANorder(mc, orderonly = T,MSTorder = c(3,2,1))
saveRDS(ord2,here('covid','data','GSE155673_pbmc','useDiffusionMap','temra','pseudotime_order.rds'))

