rm(list=ls())
library(here)
library(ggplot2)
# setwd(here())
# setwd('/Users/wenpinhou/Dropbox/')
setwd('/home/whou10/scratch16/whou10/')

source('trajectory_variability/function/01_function.R')
path = 'monocyte'

##### read in significant genes
siggene <- readRDS(paste0('trajectory_variability/hca/real/testvar/plot/perf/', path,'_venn_significant_genes.rds'))
str(siggene)
Lamian.pm = siggene[['Lamian.pm']]
Lamian.chisq = siggene[['Lamian.chisq']]
limma = siggene[['limma']]
tradeSeq = unique(unlist(siggene[['tradeSeq']]))
condiments = siggene[['condiments']]
phenopath = siggene[['phenopath']]
monocle2TrajtestCorr = siggene[['monocle2trajtestcorr']]
monocle2TrajTest = siggene[['monocle2trajtest']]

library(ggvenn)
pdf(paste0('/Users/wenpinhou/Dropbox/trajectory_variability/hca/real/testvar/plot/venn/venn_diagramm_', path, '_4methods_1.pdf'), width = 5, height = 4)
ggvenn(
  list(Lamian.pm=Lamian.pm, limma = limma, condiments=condiments, phenopath= phenopath), 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#FBB4AE"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

library(venn)
pdf(paste0('/Users/wenpinhou/Dropbox/trajectory_variability/hca/real/testvar/plot/venn/venn_diagramm_', path, '_5methods.pdf'), width = 5, height = 4)
venn(list(Lamian.pm=Lamian.pm, limma = limma, tradeSeq=tradeSeq, condiments = condiments, phenopath = phenopath), ilab=TRUE, zcolor = "style")
dev.off()



