rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
for (comparison in c('Mod_Mi', 'HD_Mi', 'Mod_Se', 'HD_Se')){
  pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/EM_pm/', comparison, '/')
  Res <- readRDS(paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/', comparison, '/cellPropTest_res.rds'))
  dir.create(pdir, recursive = T, showWarnings = F)
  names(Res)
  pdf(paste0(pdir, 'cellPropTest.pdf'), width = 2.8, height = 2.2)
  plotGene(testobj = Res, gene='prop', variable = 'type',cellProp = TRUE,  variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = F, line.alpha = 1, line.size = 0.1, point.alpha=0.5, point.size=0.5, continuous = F, sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = T, y.lab = 'Cell density')
  dev.off()
  
}



