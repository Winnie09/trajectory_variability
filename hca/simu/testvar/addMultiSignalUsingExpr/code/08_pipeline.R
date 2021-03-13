rm(list = ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
pdir <- 'demo/plot/'
rdir <- 'demo/result/'
# ##### gold standard
# trendgene <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds')
# meangene <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds')
# meantrendgene <-readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds')
# gs <- c(trendgene, meangene, meantrendgene)
## 
# ---------------------
# prepare data and test
# ---------------------
Res <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/EM_chisq/4.rds')
pdir = 'hca/simu/testvar/addMultiSignalUsingExpr/plot/EM_chisq/4/'
rdir = 'hca/simu/testvar/addMultiSignalUsingExpr/result/EM_chisq/4/'
dir.create(rdir, recursive = T)
dir.create(pdir, recursive = T)

# Res <- readRDS(paste0(rdir, 'myRes_testvar.rds'))
statistics <- Res$statistics
diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])
## if the above test works well, then refine the following codes
## --------------
## population fit
## --------------

Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')

## -----------
## clustering
## -----------
Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene)
Res$cluster <- clusterGene(Res, gene = diffgene, type = 'variable', k=3)

## --------------
## save diff gene
## --------------
allg <- diffgene
res <- data.frame(gene = allg, statistics[allg, ], cluster = Res$cluster[allg])
res <- res[order(res[, grep('^fdr.*overall$', colnames(res))]), ]
write.csv(res, paste0(rdir, 'differential_genes.csv'))

## ----------------
## plotClusterDiff
## ----------------
pdf(paste0(pdir, '/cluster_diff.pdf'), width = 3, height = 2)
plotClusterDiff(testobj=Res, gene = allg)
dev.off()

## ---------------
## plotClusterMean
## ----------------
pdf(paste0(pdir, '/cluster_mean.pdf'), width = 5, height = 3.5)
plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'variable')
dev.off()

## -----------
## GO analysis
## -----------
goRes <- GOEnrich(testobj = Res, type = 'variable')  
nn <- sapply(1:length(goRes), function(i){
  tmp <- goRes[[i]]
  tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
  # write.csv(tmp, paste0(rdir, 'cluster', i, '_GO.csv'))
  print(str(tmp))
  return(0)
})

## -----------------------
## plotClusterMeanAndDiff
## -----------------------
pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 4, height = 5.5)
print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
dev.off()

for (i in 1:max(Res$cluster)){
  print(i)
  gene <- rownames(res)[res$cluster == i]
  png(paste0(pdir, 'diffgene_groupFit_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
  print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)]))
  dev.off()
}

for (i in 1:max(Res$cluster)){
  print(i)
  gene <- rownames(res)[res$cluster == i]
  png(paste0(pdir, 'diffgene_groupDiff_cluster', i, '.png'), width = 500, height = 500, res = 200)
  print(plotClusterDiff(testobj = Res, gene = gene[1:min(length(gene), 100)], each = TRUE, sep = ':.*'))
  dev.off()
}

# --------------------------------------
# compare original and fitted expression
# --------------------------------------
png(paste0(pdir, 'fitHm.png'),width = 4000,height = 2500,res = 300)
plotFitHm(Res, type = 'variable')
dev.off()

png(paste0(pdir, 'fitHm_rownames.png'),width = 12000,height = 10000,res = 300)
print(plotFitHm(Res, type='variable', showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
dev.off()



