rm(list = ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
# ---------------------
# prepare data and test
# ---------------------
pdir = 'hca/simu/testvar/addMultiSignalUsingExpr/plot/EM_pm/4/'
Res <- readRDS(paste0(pdir, 'Res_with_clu.rds'))
statistics <- Res$statistics
diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])

## ----------------
## plotClusterDiff
## ----------------
pdf(paste0(pdir, '/cluster_diff.pdf'), width = 3, height = 2)
plotClusterDiff(testobj=Res, gene = diffgene)
dev.off()

## ---------------
## plotClusterMean
## ----------------
pdf(paste0(pdir, '/cluster_mean.pdf'), width = 5, height = 3.5)
plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'variable')
dev.off()

## -----------------------
## plotClusterMeanAndDiff
## -----------------------
pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 4, height = 5.5)
print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
dev.off()

for (i in 1:max(Res$cluster)){
  print(i)
  gene <- names(Res$cluster)[Res$cluster == i]
  png(paste0(pdir, 'diffgene_groupFit_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
  print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)]))
  dev.off()
}

for (i in 1:max(Res$cluster)){
  print(i)
  gene <- names(Res$cluster)[Res$cluster == i]
  png(paste0(pdir, 'diffgene_groupDiff_cluster', i, '.png'), width = 500, height = 500, res = 200)
  print(plotClusterDiff(testobj = Res, gene = gene[1:min(length(gene), 100)], each = TRUE, sep = ':.*'))
  dev.off()
}

# # --------------------------------------
# # compare original and fitted expression
# # --------------------------------------
# png(paste0(pdir, 'fitHm.png'),width = 4000,height = 2500,res = 300)
# plotFitHm(Res, type = 'variable')
# dev.off()
# 
# png(paste0(pdir, 'fitHm_rownames.png'),width = 12000,height = 10000,res = 300)
# print(plotFitHm(Res, type='variable', showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
# dev.off()
# 
# 


