lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/Mod_Mi/numeric_res.rds'))[[1]]
lamian <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]
limma <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/limma/Mod_Mi.rds')
limma <- rownames(limma)[limma$adj.P.Val < 0.05]
tradeseq <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/tradeSeq/Mod_Mi/testvar_res.rds')
tradeseq <- unique(unlist(sapply(tradeseq,function(i) {rownames(i)[i$adj.P.Val < 0.05]})))
tradeseq <- tradeseq[!is.na(tradeseq)]

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(tradeseq, lamian,limma),
  category.names = c("tradeSeq", "Lamian" , "Limma"),
  filename = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/perf/venn.png',
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
  rotation = 1
)


