d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/tradeSeq/monocyte/gender/sce.rds')
counts <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/build_from_tree_variability/result/monocyte/input_expr.rds')
ag <- sapply(c('BCLAF1','C19orf38','CHPT1','S100A6','TMEM69','CREBRF','DCAF6','IER3','IFNGR2'),function(i) grep(paste0('^',i,':'),rownames(counts),value=T))
library(tradeSeq)
library(RColorBrewer)
for (g in ag) {
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/tradeSeq/monocyte/gender/plot/',g,'.pdf'), width = 2.9, height = 2)
  plot(plotSmoothers(d, counts, gene = g, size = 0.1, lwd = 1))
  dev.off()  
}


