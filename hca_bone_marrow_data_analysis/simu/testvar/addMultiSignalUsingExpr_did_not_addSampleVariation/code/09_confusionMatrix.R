rm(list = ls())
library(here)
setwd(here())
library(caret)
library(pheatmap)
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
# ---------------------
# prepare data and test
# ---------------------

selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds')## trendonly
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds') ## meanonly
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds') ## trendmean

DDGType.unique <- c('trendSig', 'meanSig', 'bothSig', 'other', 'nonDDG')
signalType.unique <- c("trend only", "mean only", "trend & mean", "non-gs")


for (i in 4:1){
  print(i)
  Res <- readRDS(paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/EM_pm/', i, '.rds'))
  pdir = paste0('hca/simu/testvar/addMultiSignalUsingExpr/plot/EM_pm/', i, '/')
  dir.create(pdir, recursive = T)
  non_gs <- setdiff(rownames(Res$statistics), c(selgene1, selgene2, selgene3))
  
  statistics <- Res$statistics
  diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])
  DDGType <- getDDGType(Res)
  signalType <- c(rep('trend only', length(selgene1)), 
                  rep('mean only', length(selgene2)), 
                  rep('trend & mean', length(selgene3)), 
                  rep('non-gs', length(non_gs)))
  names(signalType) <- c(selgene1, selgene2, selgene3, non_gs)
  
  m <- confusionMatrix(data = factor(DDGType, levels = c(DDGType.unique, signalType.unique)), 
                       reference = factor(signalType, levels = c(DDGType.unique, signalType.unique)), 
                       dnn = c('DDGType','signalType'))  ## dnn = c("Prediction", "Reference")
  tb <- t(m$table[1:5, 6:9])
  tb <- tb[, DDGType.unique]
  
  tb1 <- sweep(tb, 2, colSums(tb), '/')
  pdf(paste0(pdir, '/confusionMatrix_hm_colsum1.pdf'), width = 3, height = 2.6)
  pheatmap(tb1, cluster_rows = F, cluster_cols = F)
  dev.off()

  tb2 <- tb/rowSums(tb)
  pdf(paste0(pdir, '/confusionMatrix_hm_rowsum1.pdf'), width = 3, height = 2.6)
  pheatmap(tb2, cluster_rows = F, cluster_cols = F)
  dev.off()
}

