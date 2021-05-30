rm(list=ls())
library(here)
library(ggplot2)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/EM_pm/'
pdir <- 'hca/real/testvar/plot/EM_pm/'

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')

for (path in c('erythroid', 'lymph', 'monocyte')){
  print(path)
  ## read data
  Res <- readRDS(paste0(ddir, path, '/gender/gender_res.rds'))
  res <- Res$statistics
  ## identify DDGType 
  DDGType <- getDDGType(Res)
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
  DDGType[gene]
  
  if (length(gene) > 0){
    
    pdf(paste0(pdir, path, '/gender/true_chrX.pdf'), width = 8, height = 4)
    plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
    dev.off()
  }
    
  
  print('chrY:')  
  print(diffgene.full[diffgene %in% u2])
  gene = diffgene.full[diffgene %in% u2]
  DDGType[gene]
  if (length(gene) > 0){
    pdf(paste0(pdir, path, '/gender/true_chrY.pdf'), width = 8, height = 4)
    plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
    dev.off()
  }
   
  gene = rownames(res[(nrow(res)-20): (nrow(res)), ])
  pdf(paste0(pdir, path, '/gender/nonDDG_example.pdf'), width = 8, height = 8)
  plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
  dev.off()
  if (path == 'lymph'){
    gene = c("ZFX:ENSG00000005889", "EIF1AY:ENSG00000198692", rownames(res)[grepl('PTDSS1', rownames(res))])
    pdf(paste0(pdir, path, '/gender/gene_example.pdf'), width = 5.8, height = 1.8)
    plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', palette = 'Set1', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3)
    dev.off()
  
  } else if (path == 'monocyte'){
    gene = c("XIST:ENSG00000229807",'EIF1AX:ENSG00000173674', "DDX3Y:ENSG00000067048", "EIF1AY:ENSG00000198692", rownames(res)[grepl('UQCRC1', rownames(res))], "PXN-AS1:ENSG00000255857")
    pdf(paste0(pdir, path, '/gender/gene_example.pdf'), width = 4.2, height = 5.2)
    plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', continuous = F, point.alpha = 0.1, point.size = 0.05, line.size = 0.3, palette = 'Dark2', ncol = 2)
    dev.off()
    
    gene = c('EIF1AX:ENSG00000173674', "EIF1AY:ENSG00000198692", "JPX:ENSG00000225470", "NKAPP1:ENSG00000233382", rownames(res)[grepl('UQCRC1', rownames(res))], "PXN-AS1:ENSG00000255857")
    DDGType[gene]
 #     EIF1AX:ENSG00000173674  EIF1AY:ENSG00000198692     JPX:ENSG00000225470 
 #              "bothSig"               "bothSig"               "meanSig" 
 # NKAPP1:ENSG00000233382  UQCRC1:ENSG00000010256 PXN-AS1:ENSG00000255857 
 #             "trendSig"                "nonDDG"                "nonDDG" 
    pdf(paste0(pdir, path, '/gender/gene_example_all_types.pdf'), width = 3.8, height = 4.6)
    plotGene(Res, gene = gene, variable = 'gender', plot.point = T, sep = ':.*', continuous = F, point.alpha = 0.1, point.size = 0.01, line.size = 0.3, palette = 'Dark2', ncol = 2)
    dev.off()
    
  } 
}



# 'monocyte'
# [1] "XIST:ENSG00000229807"   "JPX:ENSG00000225470"   
# [3] "RPS4X:ENSG00000198034"  "EIF1AX:ENSG00000173674"
# [5] "NKAPP1:ENSG00000233382"
#   XIST:ENSG00000229807    JPX:ENSG00000225470 
#              "bothSig"              "meanSig" 
#  RPS4X:ENSG00000198034 EIF1AX:ENSG00000173674 
#                "other"              "bothSig" 
# NKAPP1:ENSG00000233382 
#             "trendSig" 
# [1] "chrY:"
#  [1] "RPS4Y1:ENSG00000129824" "EIF1AY:ENSG00000198692"
#  [3] "TMSB4Y:ENSG00000154620" "DDX3Y:ENSG00000067048" 
#  [5] "UTY:ENSG00000183878"    "TTTY14:ENSG00000176728"
#  [7] "TTTY15:ENSG00000233864" "KDM5D:ENSG00000012817" 
#  [9] "ZFY:ENSG00000067646"    "USP9Y:ENSG00000114374" 
# RPS4Y1:ENSG00000129824 EIF1AY:ENSG00000198692 
#                "other"              "bothSig" 
# TMSB4Y:ENSG00000154620  DDX3Y:ENSG00000067048 
#              "bothSig"              "bothSig" 
#    UTY:ENSG00000183878 TTTY14:ENSG00000176728 
#              "bothSig"             "trendSig" 
# TTTY15:ENSG00000233864  KDM5D:ENSG00000012817 
#              "bothSig"              "meanSig" 
#    ZFY:ENSG00000067646  USP9Y:ENSG00000114374 
#              "bothSig"              "bothSig" 



