library(here)
setwd(here())
lamian <- readRDS(paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/Mod_Mi/numeric_res.rds'))[[1]]
lamian <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]
limma <- readRDS('covid/Su_2020_Cell/testvar/useSaverImp/result/limma/Mod_Mi.rds')
limma <- rownames(limma)[limma$adj.P.Val < 0.05]
cond <- readRDS('covid/Su_2020_Cell/testvar/useSaverImp/result/condiments/Mod_Mi/cond_gene_res.rds')
cond[,4] <- p.adjust(cond[,3])
colnames(cond)[4] <- 'FDR'
cond <- rownames(cond)[cond[,4]<0.05]

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(cond, lamian,limma),
  category.names = c("condiments", "Lamian" , "Limma"),
  filename = 'covid/Su_2020_Cell/testvar/useSaverImp/plot/perf/venn_condiments_lamian_limma.png',
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



