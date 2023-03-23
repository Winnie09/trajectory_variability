rm(list=ls())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('function/01_function.R')
pdir <- 'hca/simu/testvar/nullsimu_addSampleVariation/plot/'
Res = readRDS('hca/simu/testvar/nullsimu_addSampleVariation/result/EM_pm/res.rds')

stat = Res$statistics
str(stat)

g = rownames(stat)[grep('RBM5', rownames(stat))]
g = rownames(stat)[grep('HBB', rownames(stat))]
g = rownames(stat)[grep('TFPI', rownames(stat))]

png(paste0(pdir, sub(':.*', '', g), 'examplegene_point.png'), width = 550, height = 420, res = 200)
plotGene(Res, g, plot.point = T, sep = ':.*', line.alpha = 0, point.size = 0.05)
dev.off()

png(paste0(pdir, sub(':.*', '', g), 'examplegene_point_group.png'), width = 550, height = 420, res = 200)
plotGene(Res, g, plot.point = T, variable = 'group', sep = ':.*', line.alpha = 0, point.size = 0.05, palette = 'Dark2', continuous = F)
dev.off()


png(paste0(pdir, sub(':.*', '', g), 'examplegene_curve.png'), width = 550, height = 420, res = 200)
plotGene(Res, g, plot.point = F, variable = 'group', sep = ':.*', line.size =0.3, continuous = F)
dev.off()

pdf(paste0(pdir, sub(':.*', '', g), 'examplegene_population.pdf'), width = 3.3, height = 2.2)
plotGenePopulation(Res, g, type = 'variable', sep = ':.*', line.size = 2, palette = 'Dark2')
dev.off()

stat = stat[order(stat[,1], -stat[,3]), ]
for (g in rownames(stat)[1:30]){
  print(g)
  png(paste0(pdir, sub(':.*', '', g), 'examplegene_point.png'), width = 550, height = 420, res = 200)
  plotGene(Res, g, plot.point = T, sep = ':.*', line.alpha = 0, point.size = 0.05)
  dev.off()
  
  png(paste0(pdir, sub(':.*', '', g), 'examplegene_point_group.png'), width = 550, height = 420, res = 200)
  plotGene(Res, g, plot.point = T, variable = 'group', sep = ':.*', line.alpha = 0, point.size = 0.05, palette = 'Dark2')
  dev.off()
  
  
  png(paste0(pdir, sub(':.*', '', g), 'examplegene_curve.png'), width = 550, height = 420, res = 200)
  plotGene(Res, g, plot.point = F, variable = 'group', sep = ':.*', line.size =0.3)
  dev.off()
  
  pdf(paste0(pdir, sub(':.*', '', g), 'examplegene_population.pdf'), width = 3.2, height = 2.2)
  plotGenePopulation(Res, g, type = 'variable', sep = ':.*', line.size = 2, palette = 'Dark2')
  dev.off()
}
