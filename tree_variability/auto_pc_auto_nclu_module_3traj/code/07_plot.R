rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
library(here)
setwd(here())
suppressMessages(library(igraph))
source("/Users/wenpinhou/Dropbox/trajectory_variability/function/01_function.R")
# --------------------------------------------------------------
# input: seurat integrated object including:
# low dim reduction: umap, pca, or phate
# celltype: a dataframe, col 1 is cell name, col 2 is cell type (at least for the cells with origin cell type), col 3 is sample name
# origin: the origin cell type
# --------------------------------------------------------------
# read in data

umap = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/integrate/ser/umap.rds')
pca <- as.matrix(umap@reductions$pca@cell.embeddings)
a = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/ct/sc.rds')
ct = data.frame(cell = names(a), celltype = a, sample = sapply(names(a), function(i) sub(':.*', '', i)), stringsAsFactors = FALSE)

plotdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/auto_pc_auto_nclu_module_3traj/plot/'
rdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/auto_pc_auto_nclu_module_3traj/result/'
res <- readRDS(paste0(rdir, 'infer_tree_structure_res.rds'))

saver <- readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/matrix/saver.rds')
saver <- saver[, rownames(pca)]
str(saver)
gene <- c('CD34', 'HBB', 'CD14', 'CD3D', 'CD19', 'CD27')
gene <- sapply(gene, function(i){
  rownames(saver)[grep(paste0('^', i, ':'), rownames(saver))]
})
gene <- gene[sapply(gene, length) > 0]
plist <- list()
for (g in gene){
  expr <- saver[g, ]
  pd <- data.frame(pc1 = pca[,1], pc2 = pca[,2], expr = expr, stringsAsFactors = F)
  mycolor <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')[-c(11)]))(100)
  mycolor <- c(mycolor, rep(mycolor[100], 100))
  plist[[g]] <- 
    ggplot(data = pd, aes(x = pc1, y = pc2, color = expr)) +
    geom_scattermore() +
    scale_color_gradientn(colors = mycolor)+
    theme_classic() + ggtitle(sub(':.*', '', g)) +
    theme(legend.key.width=unit(0.2,"cm"), 
          legend.key.height=unit(0.38,"cm"), 
          axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) + labs(color = '')+
    xlab('') + ylab('')
}
png(paste0(plotdir, 'markergene.png'), res = 200, width = 1200, height = 580)
grid.arrange(grobs = plist, nrow = 2)
dev.off()

pd = data.frame(pc1 = pca[,1], pc2 = pca[,2], time = as.numeric(res$pseudotime[rownames(pca)]))
png(paste0(plotdir, 'pseudotime.png'), width = 600, height = 400, res=200)  
mycolor <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(1000)
ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) +
  geom_scattermore() +
  scale_color_gradientn(colors = mycolor)+
  theme_classic()+
  theme(legend.key.width=unit(0.2,"cm"), 
        legend.key.height=unit(0.38,"cm")) + labs(color = '')+
  xlab('Principal component 1') + ylab('Principal component 2')
dev.off()
