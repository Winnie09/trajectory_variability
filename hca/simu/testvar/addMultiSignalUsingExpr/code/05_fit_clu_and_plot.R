 rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(caret)
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
# source('/Users/wenpinhou/Dropbox/resource/myfunc/01_function.R')
source('function/01_function.R')
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- apply(data,1,sd)
  (data - cm) / csd
}

for (dataType  in seq(1,4)){
  print(dataType)
  ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/data/'
  rdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/cluster/',dataType)
  pdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/plot/', dataType)
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  Res <- readRDS(paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/EM_NOT_centered/',dataType,'.rds'))
  plotGenePopulation(Res, rownames(Res$statistics)[1], type = 'variable')
  
  fit <- getPopulationFit(Res, gene = rownames(Res$statistics)[Res$statistics[,7]<0.05], type = 'variable') ## can't run through, double check
  saveRDS(fit, paste0(rdir, '/population_fit.rds')) ##########
  Res$populationFit <- fit
  
  ## -----------
  ## clustering
  ## -----------
  Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = names(Res$fdr[Res$fdr < 0.05]))
  clu.true = readRDS(paste0(ddir, '/null/geneCluster.rds'))
  clu <- clusterGene(Res, gene = names(Res$fdr[Res$fdr < 0.05]), type = 'variable', k=max(clu.true))
  clu <- sort(clu)
  table(clu)
  saveRDS(clu, paste0(rdir, '/cluster.rds'))

  ## ---------------
  ## plotClusterDiff
  ## ----------------
  pdf(paste0(pdir, '/cluster_diff.pdf'), width = 3, height = 2)
  plotClusterDiff(testobj=Res, gene = names(Res$fdr[Res$fdr < 0.05]), cluster = clu)
  dev.off()
  
  ## ---------------
  ## plotClusterMean
  ## ----------------
  pdf(paste0(pdir, '/cluster_mean.pdf'), width = 5, height = 3.5)
  plotClusterMean(testobj=Res,cluster = clu)
  dev.off()
  
  ## ----------------
  ## confusion matrix
  ## ----------------
  ### compare signal clusters and diff clusters
  fromgene = readRDS(paste0(ddir, '/fromgene/', dataType, '.rds'))
  selgene = readRDS(paste0(ddir, '/selgene/selgene.rds'))
  clu.true = clu.true[fromgene]
  clu.true = paste0('cluster', clu.true)
  names(clu.true) <- selgene
  m <- confusionMatrix(data = as.factor(paste0('cluster',clu)), reference = as.factor(clu.true[names(clu)]), dnn = c('diffGeneCluster','SignalCluster'))  ## dnn = c("Prediction", "Reference")
  tb <- t(m$table)
  tb <- tb/rowSums(tb)
  pdf(paste0(pdir, '/confusionMatrix_hm.pdf'), width = 3.5, height = 3)
  print(pheatmap(tb))
  dev.off()
  ## compare signal type and diff type
  selgene.trend = readRDS(paste0(ddir, '/selgene/selgene1.rds'))
  selgene.mean = readRDS(paste0(ddir, '/selgene/selgene2.rds'))
  selgene.both = readRDS(paste0(ddir, '/selgene/selgene3.rds'))
  signalType <- c(rep('trendOnly', length(selgene.trend)), rep('meanOnly', length(selgene.mean)), rep('both', length(selgene.both)), rep('null', nrow(Res$statistics)-length(selgene)))
  names(signalType) <- c(selgene.trend, selgene.mean, selgene.both, setdiff(rownames(Res$statistics), selgene))
  diffType <- getDiffType(Res)
  
  tb <- sapply(sort(unique(diffType)), function(i){
    sapply(sort(unique(signalType)), function(j){
      sum(diffType == i & signalType == j)
    })
  })
  tb <- tb/rowSums(tb)
  pdf(paste0(pdir, '/confusionMatrix_hm_signalType_diffType.pdf'), width = 3.5, height = 3)
  print(pheatmap(tb, cluster_cols = FALSE, cluster_rows = FALSE))
  dev.off()
}


