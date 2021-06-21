rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
s = list()
for (comparison in c('Mod_Mi', 'HD_Mi', 'Mod_Se', 'HD_Se')){
  pdir <- paste0('covid/Su_2020_Cell/testtime/plot/EM_pm/cellPropTest/', comparison, '/')
  Res <- readRDS(paste0('covid/Su_2020_Cell/testtime/result/EM_pm/cellPropTest/', comparison, '/cellPropTest_res.rds'))
  dir.create(pdir, recursive = T, showWarnings = F)
  names(Res)
  pdf(paste0(pdir, 'cellPropTest.pdf'), width = 4.8, height = 2.2)
  plotGene(testobj = Res, gene='prop', cellProp = TRUE, free.scale = TRUE, facet.sample = FALSE, plot.point = F, line.alpha = 1, line.size = 0.1, point.alpha=0.5, point.size=0.5, continuous = T, sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = T, y.lab = 'Cell density')
  dev.off()

  pdf(paste0(pdir, 'cellPropTest_population.pdf'), width = 2.8, height = 2)
  plotGenePopulation(Res, 'prop')
  dev.off()
  
  print(comparison)
  print(Res[[1]])
  s[[comparison]] = Res[[1]]
}
s= do.call(rbind, s)
write.csv(s,'covid/Su_2020_Cell/testtime/plot/EM_pm/cellPropTest/cellPropTest_statistics.csv')


