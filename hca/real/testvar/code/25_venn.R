rm(list=ls())
library(here)
library(ggplot2)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')

lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/monocyte/gender/gender_res.rds'))[[1]]
lamian <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]
limma <- readRDS(paste0('hca/real/testvar/result/limma/monocyte/gender_res.rds'))
limma <- rownames(limma)[limma$adj.P.Val < 0.05]
tradeseq <- readRDS(paste0('hca/real/testvar/result/tradeSeq/monocyte/gender/testvar_res.rds'))
tradeseq <- unique(unlist(sapply(tradeseq,function(i) {rownames(i)[i$adj.P.Val < 0.05]})))
tradeseq <- tradeseq[!is.na(tradeseq)]

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel1")

library(VennDiagram)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/venn/venn_diagramm.pdf', width = 2.1, height = 2.1)
grid.draw(venn.diagram(
  x = list(lamian,limma,tradeseq),
  category.names = c("Lamian" , "Limma" , "tradeSeq"),
  filename = NULL,
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 2
))
dev.off()
