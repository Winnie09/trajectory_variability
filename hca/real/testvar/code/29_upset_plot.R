rm(list=ls())
library(here)
library(ggplot2)
library(RColorBrewer)
# setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
setwd('/home/whou10/scratch16/whou10/trajectory_variability/')
celltype = 'monocyte'
for (celltype in c("erythroid", "lymph", "monocyte")){
  print(celltype)
  siggene <- readRDS(paste0('hca/real/testvar/plot/perf/',celltype,'_venn_significant_genes.rds'))
  siggene[['tradeSeq']] <- unique(unlist(siggene[['tradeSeq']]))

  allg <- unique(unlist(siggene))
  mat <- sapply(names(siggene), function(met){
    print(met)
    tmp <- siggene[[met]]
    v <- rep(0, length(allg))
    names(v) <- allg
    v[tmp] <- 1
    v
  }, simplify = T)
  
  mat = do.call(cbind, mat)
  mat = as.data.frame(mat)
  library(UpSetR)
  pdf(paste0('hca/real/testvar/plot/perf/',celltype,'_num_siggene_upset_plot.pdf'),width=8.5,height=3.5)
  print(upset(mat, sets = names(mat), sets.bar.color = "#56B4E9",
        order.by = "freq", empty.intersections = "on"))
  dev.off()
}

