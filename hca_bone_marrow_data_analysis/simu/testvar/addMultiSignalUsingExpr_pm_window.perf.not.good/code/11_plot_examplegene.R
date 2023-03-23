rm(list = ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
pdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr_pm_window/plot/EM_pm/4/'
# ---------------------
# prepare data and test
# ---------------------
Res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/EM_pm/4.rds')
dir.create(pdir, recursive = T)

# Res <- readRDS(paste0(rdir, 'myRes_testvar.rds'))
statistics <- Res$statistics
diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])
## if the above test works well, then refine the following codes
## --------------
## population fit
## --------------

Res$populationFit <- getPopulationFit(Res, gene = rownames(Res$expr), type = 'variable')
v = rownames(Res$expr)
ag = c('NUCB2', 'LMNB1', 'NSMCE3', 'POLR3F')
ag = sapply(ag, function(i){
  v[grepl(i, v)]
})

pdf(paste0(pdir, 'examplegene_curve.pdf'), width = 4, height = 3.8)
plotGene(testobj = Res, gene = ag, variable = 'group', variable.text = NULL, x.lab = 'Pseudotime', y.lab = 'Expression', plot.point = F, use.palette = T, ncol = 2, sep = ':.*')
dev.off()

pdf(paste0(pdir, 'examplegene_point.pdf'), width = 4, height = 3.8)
plotGene(testobj = Res, gene = ag, variable = 'group', variable.text = NULL, x.lab = 'Pseudotime', y.lab = 'Expression', plot.point = T, use.palette = T, ncol = 2, sep = ':.*', line.alpha = 0, point.size = 0.01, point.alpha = 0.8)
dev.off()


ag = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/selgene/selgene1.rds')
png(paste0(pdir, 'selgene1.png'), width = 5000, height = 5000)
plotGene(testobj = Res, gene = ag[1:100], variable = 'group', variable.text = NULL, x.lab = 'Pseudotime', y.lab = 'Expression', plot.point = F, use.palette = T, ncol = 10, sep = ':.*')
dev.off()


ag = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/selgene/selgene2.rds')
png(paste0(pdir, 'selgene2.png'), width = 5000, height = 5000)
plotGene(testobj = Res, gene = ag[1:30], variable = 'group', variable.text = NULL, x.lab = 'Pseudotime', y.lab = 'Expression', plot.point = F, use.palette = T, ncol = 5, sep = ':.*')
dev.off()


ag = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/selgene/selgene3.rds')
png(paste0(pdir, 'selgene3.png'), width = 5000, height = 5000)
plotGene(testobj = Res, gene = ag[1:30], variable = 'group', variable.text = NULL, x.lab = 'Pseudotime', y.lab = 'Expression', plot.point = F, use.palette = T, ncol = 5, sep = ':.*')
dev.off()

