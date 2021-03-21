rm(list=ls())
library(here)
library(ggplot2)
# setwd(here())
setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/EM_pm/'
pdir <- 'hca/real/testvar/plot/EM_pm/'

## read in gold standard Sex difference genes (chrX, chrY)
# u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
# u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')

for (path in c('erythroid', 'lymph', 'monocyte')){
  print(path)
  ## read data
  Res <- readRDS(paste0(ddir, path, '/gender/gender_res.rds'))
  res <- Res$statistics
  ## identify DEGType 
  DEGType <- getDEGType(Res)
  # ------------------
  #  Evaluation Gender 
  # ------------------
  ## TP
  res <- res[order(res[,'fdr.overall'], -res[,'z.overall']), ]
  allg <- sub(':.*', '', rownames(res))
  diffgene = sub(':.*', '', rownames(res[res[,'fdr.overall'] < 0.05,]))
  diffgene.full = rownames(res[res[,'fdr.overall'] < 0.05,])
  str(diffgene)
  Res$design
  print('chrX:')
  print(diffgene.full[diffgene %in% u1])
  gene = diffgene.full[diffgene %in% u1]
  pdf(paste0(pdir, path, '/gender/true_chrX.pdf'), width = 8, height = 4)
  plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
  dev.off()
  
  print('chrY:')  
  print(diffgene.full[diffgene %in% u2])
  gene = diffgene.full[diffgene %in% u2]
  pdf(paste0(pdir, path, '/gender/true_chrY.pdf'), width = 8, height = 4)
  plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
  dev.off()
  
  gene = rownames(res[(nrow(res)-10): (nrow(res)), ])
  plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
  if (path == 'lymph'){
    gene = c("ZFX:ENSG00000005889", "EIF1AY:ENSG00000198692", rownames(res)[grepl('PTDSS1', rownames(res))])
    pdf(paste0(pdir, path, '/gender/gene_example.pdf'), width = 5.8, height = 1.8)
    plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
    dev.off()
  
  } else if (path == 'monocyte'){
    gene = c("XIST:ENSG00000229807", "ZFY:ENSG00000067646", rownames(res)[grepl('UQCRC1', rownames(res))])
    pdf(paste0(pdir, path, '/gender/gene_example.pdf'), width = 5.8, height = 1.8)
    plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
    dev.off()
  } 
}

# [1] "lymph"
# chr [1:23] "THYN1" "SCAPER" "ZFX" "GIT1" "AL592183.1" "EIF1AY" ...
# [1] "chrX:"
# [1] "ZFX:ENSG00000005889"   "RPS4X:ENSG00000198034"
# [1] "chrY:"
# [1] "EIF1AY:ENSG00000198692"
# [1] "monocyte"
# chr [1:46] "XIST" "DCAF4" "NDUFAF5" "ZFY" "RPS4Y1" "KBTBD4" ...
# [1] "chrX:"
# [1] "XIST:ENSG00000229807"    "ARHGEF9:ENSG00000131089"
# [3] "JPX:ENSG00000225470"    
# [1] "chrY:"
# [1] "ZFY:ENSG00000067646"    "RPS4Y1:ENSG00000129824"
# [3] "KDM5D:ENSG00000012817"  "DDX3Y:ENSG00000067048" 

