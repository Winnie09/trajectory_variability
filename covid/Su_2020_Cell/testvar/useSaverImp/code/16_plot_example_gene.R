rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
ddir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'
test.method = 'EM_pm'
# for (comparison in c('HD_Se', 'HD_Mi','Mod_Se','HD_Mod', 'Se_Mi', 'Recovered_Deceased')){
# for (comparison in list.files(paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method))){
for (comparison in c('Se_Mi','Mod_Mi',  'Recovered_Deceased')){
  print(comparison)
  rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method, '/', comparison, '/')
  pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', test.method, '/', comparison, '/')
  Res <- readRDS(paste0(rdir, paste0('numeric_res_with_clu.rds')))
  
  examplegene = read.csv('covid/Su_2020_Cell/testvar/useSaverImp/plot/EM_pm/sort/examplegene.csv')
  gene <- examplegene[[comparison]]
  gene <- names(sort(Res$cluster[gene]))
  ## ----------
  ## plot DEG 
  ## ----------
  id <- sort(sample(1:ncol(Res$populationFit[[1]]), ncol(Res$expr.ori)))
  Res$populationFit[[1]] <- Res$populationFit[[1]][, id]
  Res$populationFit[[2]] <- Res$populationFit[[2]][, id]
  
  png(paste0(pdir, 'examplegene_populationFit.png'), width = 1000, height = 780, res = 200)
  plotGenePopulation(Res, gene = gene, type = 'variable', facet.grid = T, axis.text.blank = T)
  dev.off()
  
  png(paste0(pdir, 'examplegene_groupDiff.png'), width = 950, height = 780, res = 200)
  plotClusterDiff(Res, gene = gene, each = T, sep = ':.*',axis.text.blank = T)
  dev.off()
  
  png(paste0(pdir, 'examplegene_sampleFit.png'), width = 1300, height = 1000, res = 200)
  plotGene(Res, gene = gene, plot.point = T, point.size = 0.05, variable = 'type', point.alpha = 0.1, line.size = 0.2)
  dev.off()
  
  for (g in gene){
    pdf(paste0(pdir, 'examplegene_populationFit_', g,'.pdf'), width = 2.4, height = 1.5)
    plotGenePopulation(Res, gene = g, type = 'variable', ncol= 4)
    dev.off()
  }
}



