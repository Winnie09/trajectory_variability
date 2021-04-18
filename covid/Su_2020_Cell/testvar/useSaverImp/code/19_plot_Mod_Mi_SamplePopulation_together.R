rm(list=ls())
library(here)
# setwd(here())
setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
test.method = 'EM_pm'
comparison = 'Mod_Mi'
print(comparison)
rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method, '/', comparison, '/')
pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', test.method, '/', comparison, '/')

Res <- readRDS(paste0(rdir, paste0('numeric_res_with_clu.rds')))
tb <- read.csv('/Users/wenpinhou/Dropbox/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/EM_pm/Mod_Mi/differential_genes_mymarkes.csv', row.names = 1)
i = tb[1:2,1]

for (c in unique(tb[,11])){
  print(c)
  i <- tb[tb[,11]==c, 1]
  png(paste0(pdir, 'cluster', c, '_SamplePopulation.png'), width = 2000, height = 2000, res = 200)
  plotGeneSampleAndPopulation(Res, gene = i, plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.1, continuous = F, palette = 'Dark2') 
  dev.off()  
}  

glist <- list(c('JAK3', 'ZEB2'),
              c('CXCR4', 'FOS'),
              c('IFI6', 'RELA'),
              c('TBX21', 'TNF'),
              c('DNMT1', 'GNLY'),
              c('BRD4', 'ISG20'),
              c('IL27RA', 'FLI1'))

for (c in 1:length(glist)){
  i = glist[[c]]
  png(paste0(pdir, 'cluster', c, '_SamplePopulation_select.png'), width = 600, height = 300, res = 200)
  plotGeneSampleAndPopulation(Res, gene = i, plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.1, continuous = F, palette = 'Dark2') 
  dev.off()  
}

