library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))

library(mclust)
for (dataType  in seq(3,4)){
  ddir <- here('hca','data','simu','testtime','addMultiSignalUsingExpr')
  rdir <- here('hca','simu','testtime', 'addMultiSignalUsingExpr', 'result','addsignal', 'EM_pm', dataType)
  pdir <- here('hca','simu','testtime', 'addMultiSignalUsingExpr', 'plot','addsignal', 'EM_pm', dataType)
  Res <- readRDS(here('hca','simu','testtime', 'addMultiSignalUsingExpr', 'result','addsignal', 'EM_pm',paste0(dataType,'.rds')))
  dir.create(rdir, recursive = T)
  
  stat <- Res$statistics
  head(stat)
  stat <- stat[order(stat[,1], -stat[,3]), ]
  ## identify XDE genes with FDR.overall < 0.05 cutoff
  diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
  str(diffgene)
  Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'time')
  clu_louvain <- clusterGene(Res, gene = diffgene, type = 'time', k=3, method = 'louvain')
  saveRDS(clu_louvain, paste0(rdir, '/cluster_louvain.rds'))
  clu_kmeans <- clusterGene(Res, gene = diffgene, type = 'time', k=5, method = 'kmeans')
  saveRDS(clu_kmeans, paste0(rdir, '/cluster_kmeans.rds'))
  
  clu_kmeans =  readRDS(paste0(rdir, '/cluster_kmeans.rds'))
  clu_gmm <- clusterGene(Res, gene = diffgene, type = 'time', k=5, method = 'GMM')
  saveRDS(clu_gmm, paste0(rdir, '/cluster_GMM.rds'))
  
  ## ----------------
  ## plot confusion matrix
  ## ----------------
  m <- caret::confusionMatrix(data = as.factor(paste0('cluster',clu_louvain)), reference = as.factor(paste0('cluster', clu_kmeans[names(clu_louvain)])),
                              dnn = c('louvain','kmeans'))  ## dnn = c("Prediction", "Reference")
  tb <- t(m$table)
  tb <- tb/rowSums(tb)
  pdf(paste0(pdir, '/confusionMatrix_hm_louvain_kmeans.pdf'), width = 3.5, height = 3)
  print(pheatmap(tb, cluster_rows = F, cluster_cols = F))
  dev.off()
  
  m <- caret::confusionMatrix(data = as.factor(paste0('cluster',clu_gmm)), reference = as.factor(paste0('cluster', clu_kmeans[names(clu_gmm)])), 
                              dnn = c('GMM','kmeans'))  ## dnn = c("Prediction", "Reference")
  tb <- t(m$table)
  tb <- tb/rowSums(tb)
  pdf(paste0(pdir, '/confusionMatrix_hm_gmm_kmeans.pdf'), width = 3.5, height = 3)
  print(pheatmap(tb, cluster_rows = F, cluster_cols = F))
  dev.off()
  ## -----------------------------------
  ## plot fitted gene expression heatmap
  ## -----------------------------------
  library(gridExtra)
  png(paste0(pdir, '/hm_kmeans.png'), width = 4100, height = 3200, res = 400)
  Res$cluster <- clu_kmeans
  plotTDEHm(
    Res,
    subsampleCell  = F,
    showCluster = TRUE,
    type = 'time',
    cellWidthTotal = 130,
    cellHeightTotal = 200
  )
  dev.off()

  png(paste0(pdir, '/hm_louvain.png'), width = 4100, height = 3200, res = 400)
  Res$cluster <- clu_louvain
  plotTDEHm(
    Res,
    subsampleCell  = F,
    showCluster = TRUE,
    type = 'time',
    cellWidthTotal = 130,
    cellHeightTotal = 200
  )
  dev.off()

  png(paste0(pdir, '/hm_gmm.png'), width = 4100, height = 3200, res = 400)
  Res$cluster <- clu_gmm
  plotTDEHm(
    Res,
    subsampleCell  = F,
    showCluster = TRUE,
    type = 'time',
    cellWidthTotal = 130,
    cellHeightTotal = 200
  )
  dev.off()
}



