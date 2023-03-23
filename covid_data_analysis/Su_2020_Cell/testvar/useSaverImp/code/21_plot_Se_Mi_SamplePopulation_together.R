rm(list=ls())
library(here)
# setwd(here())
setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
test.method = 'EM_pm'
comparison = 'Se_Mi'
print(comparison)
rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method, '/', comparison, '/')
pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', test.method, '/', comparison, '/')

Res <- readRDS(paste0(rdir, paste0('numeric_res_with_clu.rds')))
tb <- read.csv(paste0('/Users/wenpinhou/Dropbox/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/EM_pm/', comparison, '/differential_genes_zeyu.csv'), row.names = 1)
i = tb[1:2,1]
str(tb)

for (c in unique(tb[,11])){
  print(c)
  i <- tb[tb[,11]==c, 1]
  png(paste0(pdir, 'cluster', c, '_SamplePopulation.png'), width = 2000, height = 2000, res = 200)
  plotGeneSampleAndPopulation(Res, gene = i, plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.1, continuous = F, palette = 'Dark2') 
  dev.off()  
}  

glist <- list(c('IL2RG', 'HLA-DPA1'),
              c('RELB', 'STAT3'),
              c('JUNB', 'FOS'),
              c('EOMES', 'TBX21'),
              c('ADAR', 'IKZF3'), 
              c('CD37', 'IFNGR2'),
              c('IFI35', 'S100A8'),
              c('IL10RA', 'STAT5A'))
for (c in 1:length(glist)){
  i = glist[[c]]
  png(paste0(pdir, 'cluster', c, '_SamplePopulation_select.png'), width = 600, height = 300, res = 200)
  plotGeneSampleAndPopulation(Res, gene = i, plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.1, line.alpha = 0.6, continuous = F, palette = 'Dark2') 
  dev.off()  
}

g <- c('TCF7', 'GZMB', 'IL7RA', 'KLRD1', 'PDCD1')
for (c in g){
  if (c %in% rownames(Res$covariateGroupDiff)){
    png(paste0(pdir, c, '_SamplePopulation.png'), width = 600, height = 450, res = 200)
    plotGeneSampleAndPopulation(Res, gene = c, plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.05, continuous = F, palette = 'Dark2') 
    dev.off()  
  }
}  

