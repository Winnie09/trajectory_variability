
rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
pdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/hca/simu/testvar/nullsimu/plot/'
Res = readRDS('hca/simu/testvar/nullsimu/result/EM_pm/res.rds')

stat = Res$statistics
str(stat)

g = rownames(stat)[grep('RBM5', rownames(stat))]

png(paste0(pdir, 'examplegene_point.png'), width = 550, height = 420, res = 200)
plotGene(Res, g, plot.point = T, sep = ':.*', line.alpha = 0, point.size = 0.05)
dev.off()

png(paste0(pdir, 'examplegene_point_group.png'), width = 550, height = 420, res = 200)
plotGene(Res, g, plot.point = T, variable = 'group', sep = ':.*', line.alpha = 0, point.size = 0.05, palette = 'Dark2')
dev.off()


png(paste0(pdir, 'examplegene_curve.png'), width = 550, height = 420, res = 200)
plotGene(Res, g, plot.point = F, variable = 'group', sep = ':.*', line.size =0.3)  
dev.off()

pdf(paste0(pdir, 'examplegene_population.pdf'), width = 3.2, height = 2.2)
plotGenePopulation(Res, g, type = 'variable', sep = ':.*', line.size = 2, palette = 'Dark2')
dev.off()




